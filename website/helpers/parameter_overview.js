(() => {
  const widget = document.currentScript.previousElementSibling;
  const parameters = JSON.parse(widget.querySelector('.parameter-overview-data').textContent).parameters;
  const search = widget.querySelector('input[type="search"]');
  const topic = widget.querySelector('select');
  const filters = [...widget.querySelectorAll('[data-status]')];
  const list = widget.querySelector('.parameter-topic-list');
  const summary = widget.querySelector('.parameter-overview-summary');
  const empty = widget.querySelector('.parameter-empty-state');
  let status = 'all';
  const expanded = new Set();
  const escape = value => String(value ?? '').replaceAll('&', '&amp;').replaceAll('<', '&lt;')
    .replaceAll('>', '&gt;').replaceAll('"', '&quot;').replaceAll("'", '&#39;');
  const statusLabel = value => value === 'Must specify' ? 'Required' : value;
  const statusClass = value => value === 'Must specify' ? 'required' : value.toLowerCase();
  const displayTopic = value => value === 'required' ? 'Core configuration' : value;
  const searchText = parameter => Object.values(parameter).join(' ').toLowerCase();
  const topics = [...new Set(parameters.map(parameter => parameter.topic))];
  topics.forEach(value => topic.add(new Option(displayTopic(value), value)));

  function renderRow(parameter) {
    const defaultValue = parameter.status === 'Must specify' ? 'Supply a value' : parameter.default_value;
    const detail = (label, value, code = true) => value ? `<div class="parameter-detail"><dt>${escape(label)}</dt><dd>${code ? `<code>${escape(value)}</code>` : escape(value)}</dd></div>` : '';
    return `<details class="parameter-row" data-param="${escape(parameter.param_name)}"${expanded.has(parameter.param_name) ? ' open' : ''}>
      <summary class="parameter-head">
        <span class="parameter-name-wrap"><code class="parameter-name">${escape(parameter.param_name)}</code>
          <span class="parameter-status ${statusClass(parameter.status)}">${statusLabel(parameter.status)}</span></span>
        <span class="parameter-default"><span>Default</span><code>${escape(defaultValue)}</code></span>
        <span class="parameter-description">${escape(parameter.description)}</span>
      </summary>
      <dl class="parameter-details">
        ${detail('Type', `${parameter.data_type} · ${parameter.cardinality}`, false)}
        ${detail('Allowed values', parameter.allowed_values || 'Any value matching the type', Boolean(parameter.allowed_values))}
        ${detail('Example', parameter.examples)}
        ${detail('Used by', parameter.part_of, false)}
      </dl>
    </details>`;
  }

  function render() {
    const words = search.value.trim().toLowerCase().split(/\s+/).filter(Boolean);
    const matched = parameters.filter(parameter => (!topic.value || parameter.topic === topic.value)
      && words.every(word => searchText(parameter).includes(word)));
    const visible = matched.filter(parameter => status === 'all' || parameter.status === status);
    filters.forEach(button => {
      const active = button.dataset.status === status;
      button.setAttribute('aria-pressed', String(active));
      button.classList.toggle('is-active', active);
      const count = matched.filter(parameter => button.dataset.status === 'all' || parameter.status === button.dataset.status).length;
      button.querySelector('span').textContent = count;
    });
    list.innerHTML = topics.map(value => {
      const rows = visible.filter(parameter => parameter.topic === value);
      return rows.length ? `<section class="parameter-topic"><h4>${escape(displayTopic(value))}<span>${rows.length}</span></h4>
        <div class="parameter-list">${rows.map(renderRow).join('')}</div></section>` : '';
    }).join('');
    summary.textContent = `${visible.length} of ${parameters.length} parameters`;
    empty.hidden = visible.length !== 0;
    list.querySelectorAll('details').forEach(row => row.addEventListener('toggle', () => {
      if (row.open) expanded.add(row.dataset.param); else expanded.delete(row.dataset.param);
    }));
  }
  search.addEventListener('input', render);
  topic.addEventListener('change', render);
  filters.forEach(button => button.addEventListener('click', () => { status = button.dataset.status; render(); }));
  widget.querySelector('.parameter-reset').addEventListener('click', () => {
    search.value = ''; topic.value = ''; status = 'all'; render(); search.focus();
  });
  render();
})();
