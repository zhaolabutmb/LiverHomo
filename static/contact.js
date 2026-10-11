document.addEventListener('DOMContentLoaded', () => {
  const form = document.getElementById('contact-form');
  const status = document.getElementById('contact-status');
  const submit = document.getElementById('contact-submit');
  if (!form || !status || !submit) return;

  form.addEventListener('submit', async (event) => {
    event.preventDefault();
    submit.disabled = true;
    status.textContent = 'Sending message...';
    try {
      const response = await fetch('/api/contact', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          name: document.getElementById('contact-name').value,
          email: document.getElementById('contact-email').value,
          message: document.getElementById('contact-message').value,
        }),
      });
      const payload = await response.json();
      if (!response.ok) throw new Error(payload.error || 'Message could not be sent.');
      form.reset();
      status.textContent = payload.message;
    } catch (error) {
      status.textContent = error.message || 'Message could not be sent.';
    } finally {
      submit.disabled = false;
    }
  });
});