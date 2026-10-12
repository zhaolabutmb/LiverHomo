"""Run the assistant against fixed questions and score the answers.

Usage (from dev/)::

    python scripts/eval_assistant.py --label baseline [--max-tokens 400000] [--only id1,id2]

Each run calls the OpenAI API. Token usage is printed per question and the
run stops early once --max-tokens is reached. Results are written to
scripts/eval_results/<label>.json (ignored by git).
"""
import argparse
import inspect
import json
import re
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import app as webapp  # noqa: E402

CASES_PATH = Path(__file__).with_name("assistant_eval_cases.json")
RESULTS_DIR = Path(__file__).with_name("eval_results")


class UsageRecorder:
    """Wraps urlopen so every OpenAI response's token usage and tool calls are recorded."""

    def __init__(self):
        self.calls = []
        self._urlopen = webapp.urllib_request.urlopen

    def install(self):
        recorder = self

        class Response:
            def __init__(self, inner):
                self.inner = inner

            def __enter__(self):
                return self

            def __exit__(self, *exc):
                self.inner.close()

            def read(self):
                body = self.inner.read()
                try:
                    package = json.loads(body.decode("utf-8"))
                    message = package["choices"][0]["message"]
                    recorder.calls.append({
                        "usage": package.get("usage", {}),
                        "tool_calls": [
                            {"name": call["function"]["name"], "arguments": call["function"]["arguments"]}
                            for call in message.get("tool_calls") or []
                        ],
                        "finish_reason": package["choices"][0].get("finish_reason"),
                    })
                except (ValueError, KeyError, IndexError):
                    pass
                return body

        def urlopen(*args, **kwargs):
            return Response(recorder._urlopen(*args, **kwargs))

        webapp.urllib_request.urlopen = urlopen

    def take(self):
        calls, self.calls = self.calls, []
        return calls


def score(case, answer, visualization):
    text = (answer or "").lower()
    failures = []
    for needle in case.get("expect_all", []):
        if needle.lower() not in text:
            failures.append(f"missing '{needle}'")
    if case.get("expect_any") and not any(needle.lower() in text for needle in case["expect_any"]):
        failures.append(f"none of {case['expect_any']}")
    for needle in case.get("forbid", []):
        if needle.lower() in text:
            failures.append(f"contains forbidden '{needle}'")
    expected_viz = case.get("expect_viz")
    if expected_viz:
        if not visualization:
            failures.append("no chart")
        else:
            title = visualization.get("title", "")
            if "chart_type" in expected_viz and visualization.get("chart_type") != expected_viz["chart_type"]:
                failures.append(f"chart_type {visualization.get('chart_type')}")
            if "title_startswith" in expected_viz and not title.startswith(expected_viz["title_startswith"]):
                failures.append(f"chart title '{title}'")
            if "title_contains" in expected_viz and expected_viz["title_contains"] not in title:
                failures.append(f"chart title '{title}'")
    if re.fullmatch(r"\s*i do not know\.?\s*", text):
        failures.append("answered 'I do not know'")
    return failures


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", required=True)
    parser.add_argument("--max-tokens", type=int, default=400_000)
    parser.add_argument("--only", default="")
    args = parser.parse_args()

    cases = json.loads(CASES_PATH.read_text())
    if args.only:
        wanted = set(args.only.split(","))
        cases = [case for case in cases if case["id"] in wanted]

    supports_history = "history" in inspect.signature(webapp.call_openai).parameters
    recorder = UsageRecorder()
    recorder.install()
    results, total_tokens = [], 0
    for case in cases:
        if total_tokens >= args.max_tokens:
            print(f"Token budget reached ({total_tokens}); skipping remaining cases.")
            break
        started = time.time()
        kwargs = {"history": case["history"]} if supports_history and case.get("history") else {}
        try:
            payload, status = webapp.call_openai(case["question"], **kwargs)
        except Exception as exc:  # report and continue
            payload, status = {"error": {"message": repr(exc)}}, 0
        calls = recorder.take()
        tokens = sum(call["usage"].get("total_tokens", 0) for call in calls)
        total_tokens += tokens
        answer = payload.get("answer") if status == 200 else None
        visualization = payload.get("visualization") if status == 200 else None
        if visualization:
            visualization = {key: value for key, value in visualization.items() if key != "image_base64"}
        failures = score(case, answer, visualization) if status == 200 else [f"HTTP {status}: {payload.get('error', {}).get('message', '')[:120]}"]
        results.append({
            "id": case["id"],
            "question": case["question"],
            "passed": not failures,
            "failures": failures,
            "answer": answer,
            "visualization": visualization,
            "status": status,
            "tool_calls": [tool for call in calls for tool in call["tool_calls"]],
            "finish_reasons": [call["finish_reason"] for call in calls],
            "tokens": tokens,
            "prompt_tokens": sum(call["usage"].get("prompt_tokens", 0) for call in calls),
            "completion_tokens": sum(call["usage"].get("completion_tokens", 0) for call in calls),
            "seconds": round(time.time() - started, 1),
        })
        if status not in (200, 0) and len(results) >= 2 and all(r["status"] == status for r in results[-2:]):
            print(f"Two consecutive HTTP {status} errors; stopping.")
            break
        mark = "PASS" if not failures else "FAIL"
        print(f"{mark} {case['id']:<22} {tokens:>7} tok {results[-1]['seconds']:>6}s  {'; '.join(failures)}", flush=True)

    passed = sum(result["passed"] for result in results)
    print(f"\n{args.label}: {passed}/{len(results)} passed, {total_tokens} tokens "
          f"(prompt {sum(r['prompt_tokens'] for r in results)}, completion {sum(r['completion_tokens'] for r in results)})")
    RESULTS_DIR.mkdir(exist_ok=True)
    out = RESULTS_DIR / f"{args.label}.json"
    out.write_text(json.dumps({"label": args.label, "passed": passed, "total": len(results),
                               "tokens": total_tokens, "results": results}, indent=2))
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
