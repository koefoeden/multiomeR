(() => {
  const widget = document.currentScript.previousElementSibling;
  const parameters = JSON.parse(widget.querySelector('.parameter-overview-data').textContent).parameters;
  const search = widget.querySelector('input[type="search"]');
  const tabs = widget.querySelector('[role="tablist"]');
  const panel = widget.querySelector('[role="tabpanel"]');
  const list = widget.querySelector('.parameter-section-list');
  const summary = widget.querySelector('.parameter-overview-summary');
  const empty = widget.querySelector('.parameter-empty-state');
  const workflows = [
    ['aggregation', 'Main workflow'],
    ['differential_analyses', 'Differential analyses'],
    ['genetic_enrichment', 'Genetic enrichment'],
    ['peak_gene_correlation', 'Peak–gene correlation']
  ].filter(([scope]) => parameters.some(parameter => parameter.scope === scope));
  let workflow = workflows[0][0];
  const expanded = new Set();
  const escape = value => String(value ?? '').replaceAll('&', '&amp;').replaceAll('<', '&lt;')
    .replaceAll('>', '&gt;').replaceAll('"', '&quot;').replaceAll("'", '&#39;');
  const searchText = parameter => Object.values(parameter).join(' ').toLowerCase();
  tabs.innerHTML = workflows.map(([scope, label]) => `<button type="button" role="tab" class="parameter-workflow" id="workflow-${scope}" data-scope="${scope}" aria-controls="parameter-panel">${label}<span>${parameters.filter(parameter => parameter.scope === scope).length}</span></button>`).join('');
  const buttons = [...tabs.querySelectorAll('button')];

  function renderRow(parameter) {
    const defaultValue = parameter.status === 'Required' ? 'Supply a value' : parameter.default_value;
    const detail = (label, value, code = true) => value ? `<div class="parameter-detail"><dt>${escape(label)}</dt><dd>${code ? `<code>${escape(value)}</code>` : escape(value)}</dd></div>` : '';
    return `<details class="parameter-row" id="${escape(parameter.param_name)}" tabindex="-1" data-param="${escape(parameter.param_name)}"${expanded.has(parameter.param_name) ? ' open' : ''}>
      <summary class="parameter-head">
        <span class="parameter-name-wrap">
          <span class="parameter-title"><strong>${escape(parameter.short_name)}</strong><span class="parameter-status ${parameter.status.toLowerCase()}">${parameter.status}</span></span>
          <code class="parameter-name">${escape(parameter.param_name)}</code>
        </span>
        <span class="parameter-default"><span>Default</span><code>${escape(defaultValue)}</code></span>
        <span class="parameter-description">${escape(parameter.description)}</span>
      </summary>
      <p><a href="#${encodeURIComponent(parameter.param_name)}" class="parameter-permalink">Link to this parameter</a></p>
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
    const available = parameters.filter(parameter => parameter.scope === workflow);
    const visible = available.filter(parameter => words.every(word => searchText(parameter).includes(word)))
      .sort((a, b) => a.short_name.localeCompare(b.short_name, 'en', {sensitivity: 'base'}) || a.param_name.localeCompare(b.param_name));
    buttons.forEach(button => {
      const active = button.dataset.scope === workflow;
      button.setAttribute('aria-selected', String(active));
      button.tabIndex = active ? 0 : -1;
    });
    panel.setAttribute('aria-labelledby', `workflow-${workflow}`);
    list.innerHTML = ['Required', 'Defaulted', 'Optional'].map(status => {
      const rows = visible.filter(parameter => parameter.status === status);
      return rows.length ? `<section class="parameter-section" data-status="${status}"><h2>${status}<span>${rows.length}</span></h2>
        <div class="parameter-list">${rows.map(renderRow).join('')}</div></section>` : '';
    }).join('');
    summary.textContent = `${workflows.find(([scope]) => scope === workflow)[1]}: ${visible.length} of ${available.length} parameters`;
    empty.hidden = visible.length !== 0;
    list.querySelectorAll('details').forEach(row => row.addEventListener('toggle', () => {
      if (row.open) expanded.add(row.dataset.param); else expanded.delete(row.dataset.param);
    }));
  }

  function navigate() {
    let fragment;
    try { fragment = decodeURIComponent(window.location.hash.slice(1)); } catch { return; }
    const parameter = parameters.find(parameter => parameter.param_name === fragment);
    const scope = fragment.startsWith('workflow=') ? fragment.slice('workflow='.length) : null;
    if (parameter) workflow = parameter.scope;
    else if (workflows.some(([value]) => value === scope)) workflow = scope;
    else if (!fragment) workflow = workflows[0][0];
    else return;
    search.value = '';
    if (parameter) expanded.add(parameter.param_name);
    render();
    if (parameter) {
      const row = document.getElementById(parameter.param_name);
      row.focus({preventScroll: true});
      row.scrollIntoView({block: 'center'});
    }
  }
  buttons.forEach(button => button.addEventListener('click', () => {
    const hash = `#workflow=${button.dataset.scope}`;
    if (window.location.hash === hash) navigate(); else window.location.hash = hash;
  }));
  tabs.addEventListener('keydown', event => {
    const index = buttons.indexOf(event.target);
    if (index < 0) return;
    const next = {ArrowRight: (index + 1) % buttons.length, ArrowLeft: (index + buttons.length - 1) % buttons.length, Home: 0, End: buttons.length - 1}[event.key];
    if (next === undefined) return;
    event.preventDefault();
    buttons[next].focus();
    buttons[next].click();
  });
  search.addEventListener('input', render);
  widget.querySelector('.parameter-reset').addEventListener('click', () => {
    search.value = ''; render(); search.focus();
  });
  window.addEventListener('hashchange', navigate);
  widget.addEventListener('click', event => {
    if (event.target.closest('.parameter-permalink')?.hash === window.location.hash) navigate();
  });
  render();
  navigate();
})();
