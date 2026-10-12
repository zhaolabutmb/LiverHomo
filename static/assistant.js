document.addEventListener('DOMContentLoaded', () => {
  const form = document.getElementById('assistant-form');
  const input = document.getElementById('assistant-input');
  const messages = document.getElementById('assistant-messages');
  const submit = document.getElementById('assistant-submit');
  if (!form || !input || !messages || !submit) return;
  const history = [];
  const MAX_HISTORY = 6;

  const escapeHtml = (text) => text.replace(/[&<>"']/g, (ch) => ({
    '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;',
  })[ch]);

  // Escape first, then allow only **bold** and line breaks from the model's text.
  const formatAnswer = (text) => escapeHtml(text)
    .replace(/\*\*(.+?)\*\*/g, '<strong>$1</strong>')
    .replace(/\n/g, '<br>');

  const addMessage = (role, text) => {
    const message = document.createElement('p');
    message.className = `assistant-message ${role}`;
    if (role === 'assistant') {
      message.innerHTML = formatAnswer(text);
    } else {
      message.textContent = text;
    }
    messages.appendChild(message);
    messages.scrollTop = messages.scrollHeight;
  };

  const addVisualization = (visualization) => {
    if (!visualization || !visualization.image_base64) return;
    const figure = document.createElement('figure');
    figure.className = 'assistant-visualization';
    const img = document.createElement('img');
    img.src = `data:image/png;base64,${visualization.image_base64}`;
    img.alt = visualization.title || 'Generated chart';
    figure.appendChild(img);
    if (visualization.title) {
      const caption = document.createElement('figcaption');
      caption.textContent = visualization.title;
      figure.appendChild(caption);
    }
    if (visualization.destination_url && visualization.destination_label) {
      const link = document.createElement('a');
      link.href = visualization.destination_url;
      link.textContent = visualization.destination_label;
      link.className = 'assistant-visualization-link';
      figure.appendChild(link);
    }
    messages.appendChild(figure);
    messages.scrollTop = messages.scrollHeight;
  };

  addMessage('assistant', 'Ask about genes in liver diseases, for example: "Is SPP1 upregulated in HCC macrophages?", "Top 10 genes down in PBC cholangiocytes", or "Show ALB expression across cell types".');

  form.addEventListener('submit', async (event) => {
    event.preventDefault();
    const query = input.value.trim();
    if (!query) return;

    addMessage('user', query);
    input.value = '';
    input.disabled = true;
    submit.disabled = true;
    try {
      const response = await fetch('/api/query', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ query, history: history.slice(-MAX_HISTORY) }),
      });
      const payload = await response.json();
      if (response.ok && payload.answer) {
        addMessage('assistant', payload.answer);
        addVisualization(payload.visualization);
        history.push({ role: 'user', content: query }, { role: 'assistant', content: payload.answer });
      } else {
        addMessage('assistant', payload?.error?.message || 'Sorry, something went wrong. Please try again.');
      }
    } catch (error) {
      addMessage('assistant', 'The assistant could not be reached. Please try again.');
    } finally {
      input.disabled = false;
      submit.disabled = false;
      input.focus();
    }
  });
});