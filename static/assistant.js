document.addEventListener('DOMContentLoaded', () => {
  const form = document.getElementById('assistant-form');
  const input = document.getElementById('assistant-input');
  const messages = document.getElementById('assistant-messages');
  const submit = document.getElementById('assistant-submit');
  if (!form || !input || !messages || !submit) return;

  const addMessage = (role, text) => {
    const message = document.createElement('p');
    message.className = `assistant-message ${role}`;
    message.textContent = text;
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

  addMessage('assistant', 'Ask a question about the data available in this portal.');

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
        body: JSON.stringify({ query }),
      });
      const payload = await response.json();
      addMessage('assistant', response.ok ? payload.answer || 'I do not know' : 'I do not know');
      if (response.ok) addVisualization(payload.visualization);
    } catch (error) {
      addMessage('assistant', 'I do not know');
    } finally {
      input.disabled = false;
      submit.disabled = false;
      input.focus();
    }
  });
});