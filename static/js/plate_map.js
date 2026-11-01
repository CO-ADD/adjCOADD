(function () {
  'use strict';

  const dataEl = document.getElementById('plate-map-data');
  if (!dataEl) { return; }

  let data;
  try {
    data = JSON.parse(dataEl.textContent);
  } catch (e) {
    console.error('[plate_map] Invalid plate data', e);
    return;
  }

  const gridEl = document.getElementById('plate-map-grid');
  const metricsEl = document.getElementById('plate-map-metrics');
  const legendEl = document.getElementById('plate-map-legend');
  const typeLegendEl = document.getElementById('plate-map-type-legend');
  const statusEl = document.getElementById('plate-map-status');
  const tooltipEl = document.getElementById('plate-map-tooltip');

  if (!gridEl) { return; }

  // Well type -> border colour and centre icon (single vs multiple compounds)
  const TYPE_META = {
    NEG:  { label: 'NEG',   single: 'bi-dash-circle',  multi: 'bi-dash-square',  border: '#6a0000', fill: '#8b0000', text: '#ffffff' },
    POS:  { label: 'POS',   single: 'bi-plus-circle',  multi: 'bi-plus-square',  border: '#004a00', fill: '#006400', text: '#ffffff' },
    CTRL: { label: 'CTRL',  single: 'bi-record-circle', multi: 'bi-record-btn',  border: '#8b0000', fill: '#d3d3d3', text: '#212529' },
    SMP:  { label: 'SMP',   single: 'bi-circle',       multi: 'bi-square',       border: '#808080', fill: '#d3d3d3', text: '#212529' },
    SKIP: { label: 'SKIP',  single: 'bi-slash-circle', multi: 'bi-slash-square', border: '#999999', fill: '#f2f2f2', text: '#212529' },
    EMPTY:{ label: 'Empty', single: null, multi: null, symbol: '',  border: '#000000', fill: '#ffffff', text: '#212529' },
    NA:   { label: '\u2014', single: null, multi: null, symbol: '', border: '#ced4da', fill: '#ffffff', text: '#212529' }
  };

  const FALLBACK_TYPE = (data.show_type === false) ? 'NA' : 'EMPTY';

  function rawType(w) {
    if (!w) { return FALLBACK_TYPE; }
    if (w.is_negcontrol) { return 'NEG'; }
    if (w.is_poscontrol) { return 'POS'; }
    if (w.is_control) { return 'CTRL'; }
    if (w.is_skip) { return 'SKIP'; }
    if (w.is_sample) { return 'SMP'; }
    return FALLBACK_TYPE;
  }

  // A sample well with no compound batch is rendered as Empty.
  function effectiveType(w) {
    const t = rawType(w);
    if (t === 'SMP') {
      const cmpCount = (w && w.cmpbatch) ? w.cmpbatch.length : 0;
      if (cmpCount === 0) { return 'EMPTY'; }
    }
    return t;
  }

  // Centre icon: single compound -> circle variant, multiple -> square variant.
  // Empty/NA wells show nothing.
  function typeIconHTML(w) {
    const t = effectiveType(w);
    const tm = TYPE_META[t] || TYPE_META.NA;
    const cmpCount = (w && w.cmpbatch) ? w.cmpbatch.length : 0;
    const icon = (cmpCount > 1) ? tm.multi : tm.single;
    if (icon) { return '<i class="bi ' + icon + '"></i>'; }
    return tm.symbol || '';
  }

  const metricNames = Object.keys(data.metrics || {});
  const wells = data.wells || [];

  if (wells.length === 0) {
    if (statusEl) { statusEl.textContent = 'No well data available for this plate.'; }
    return;
  }

  if (metricNames.length === 0) {
    if (statusEl) { statusEl.textContent = 'No heatmap properties available.'; }
    return;
  }

  let activeMetric = metricNames[0];

  const wellById = {};
  wells.forEach(function (w) { wellById[w.id] = w; });

  // ---------- helpers ----------
  function metricValue(w, m) {
    if (!w) { return null; }
    if (m.key === 'inhibition') {
      const v = w.inhibition;
      return (v === null || v === undefined || !isFinite(Number(v))) ? null : Number(v);
    }
    if (m.key === 'readout') {
      const arr = w.readout || [];
      return arr.length ? Number(arr[0]) : null;
    }
    if (m.key === 'conc') {
      const arr = w.conc || [];
      const idx = m.index || 0;
      return (arr.length > idx) ? Number(arr[idx]) : null;
    }
    return null;
  }

  function fmt(v) {
    if (v === null || v === undefined || !isFinite(Number(v))) { return '—'; }
    const n = Number(v);
    if (Math.abs(n) >= 100) { return n.toFixed(0); }
    if (Math.abs(n) >= 10) { return n.toFixed(1); }
    if (Math.abs(n) >= 1) { return n.toFixed(2); }
    return n.toFixed(3);
  }

  function computeStats(m) {
    const vals = [];
    wells.forEach(function (w) {
      const v = metricValue(w, m);
      if (v !== null) { vals.push(v); }
    });
    if (vals.length === 0) { return { min: 0, max: 1 }; }
    return { min: Math.min.apply(null, vals), max: Math.max.apply(null, vals) };
  }

  function colorFor(w, stats, m) {
    const v = metricValue(w, m);
    if (v === null) { return null; }

    if (m.mode === 'diverging') {
      const lo = (m.min !== undefined) ? Number(m.min) : stats.min;
      const hi = (m.max !== undefined) ? Number(m.max) : stats.max;
      const t = hi === lo ? 0.5 : Math.min(1, Math.max(0, (v - lo) / (hi - lo)));
      const hue = 120 * (1 - t); // 0 (low) -> green, 1 (high) -> red
      return 'hsl(' + hue + ', 65%, 45%)';
    }

    // sequential: light -> dark red, normalized across non-null values
    const t = (stats.max === stats.min) ? 0 : Math.min(1, Math.max(0, (v - stats.min) / (stats.max - stats.min)));
    const light = 92 - t * 57; // 92% -> 35%
    return 'hsl(0, 65%, ' + light + '%)';
  }

  function textColorFor(bg) {
    // pick readable text based on background lightness
    const m = bg.match(/hsl\([\d.]+,?\s*[\d.]+%,\s*([\d.]+)%\)/);
    if (m) { return Number(m[1]) > 65 ? '#212529' : '#ffffff'; }
    return '#212529';
  }

  function wellFlags(w) {
    switch (rawType(w)) {
      case 'NEG': return ['NEG'];
      case 'POS': return ['POS'];
      case 'CTRL': return ['CTRL'];
      case 'SMP': return ['SMP'];
      case 'SKIP': return ['SKIP'];
      default: return [];
    }
  }

  function buildTooltip(w) {
    const parts = [];
    parts.push('<div class="plate-map-tt-head">' + w.id + '</div>');
    const flags = wellFlags(w);
    if (flags.length) {
      parts.push('<div class="plate-map-tt-flags">' + flags.join(' &middot; ') + '</div>');
    }

    const cmpbatch = w.cmpbatch || [];
    const conc = w.conc || [];
    const concUnit = w.conc_unit || [];
    if (cmpbatch.length) {
      const items = cmpbatch.map(function (cb, i) {
        let txt = cb;
        if (conc[i] !== undefined && conc[i] !== null) {
          txt += ' (' + fmt(conc[i]) + (concUnit[i] ? ' ' + concUnit[i] : '') + ')';
        }
        return txt;
      });
      parts.push('<div class="plate-map-tt-row"><span>Compound(s)</span><b>' + items.join('<br>') + '</b></div>');
    } else {
      parts.push('<div class="plate-map-tt-row"><span>Compound(s)</span><b>—</b></div>');
    }

    if (w.barcode !== undefined) {
      parts.push('<div class="plate-map-tt-row"><span>Barcode</span><b>' + (w.barcode || '—') + '</b></div>');
    }
    if (w.readout !== undefined && w.readout.length) {
      parts.push('<div class="plate-map-tt-row"><span>Readout</span><b>' + w.readout.map(fmt).join(', ') + '</b></div>');
    }
    if (w.inhibition !== undefined && w.inhibition !== null) {
      parts.push('<div class="plate-map-tt-row"><span>Inhibition</span><b>' + fmt(w.inhibition) + '%</b></div>');
    }
    if (w.act_type) {
      parts.push('<div class="plate-map-tt-row"><span>Activity</span><b>' + w.act_type + '</b></div>');
    }

    return parts.join('');
  }

  function showTooltip(w, cellEl) {
    tooltipEl.innerHTML = buildTooltip(w);
    tooltipEl.hidden = false;
    const rect = cellEl.getBoundingClientRect();
    let left = rect.right + 8;
    let top = rect.top;
    const ttWidth = tooltipEl.offsetWidth;
    const ttHeight = tooltipEl.offsetHeight;
    if (left + ttWidth > window.innerWidth - 8) {
      left = rect.left - ttWidth - 8;
    }
    if (top + ttHeight > window.innerHeight - 8) {
      top = window.innerHeight - ttHeight - 8;
    }
    tooltipEl.style.left = left + 'px';
    tooltipEl.style.top = top + 'px';
  }

  function hideTooltip() {
    tooltipEl.hidden = true;
  }

  // ---------- build grid ----------
  function buildGrid() {
    gridEl.style.gridTemplateColumns = '26px repeat(' + data.n_cols + ', minmax(20px, 1fr))';
    gridEl.style.minWidth = (26 + data.n_cols * 24) + 'px';

    // top-left corner + column numbers
    const corner = document.createElement('div');
    corner.className = 'plate-map-label plate-map-corner';
    gridEl.appendChild(corner);
    for (let c = 1; c <= data.n_cols; c++) {
      const label = document.createElement('div');
      label.className = 'plate-map-label';
      label.textContent = c;
      gridEl.appendChild(label);
    }

    // rows
    for (let r = 0; r < data.n_rows; r++) {
      const rowLabel = document.createElement('div');
      rowLabel.className = 'plate-map-label';
      rowLabel.textContent = data.row_labels[r] || '';
      gridEl.appendChild(rowLabel);

      for (let c = 1; c <= data.n_cols; c++) {
        const cell = document.createElement('div');
        cell.className = 'plate-map-cell';
        cell.dataset.wellId = '';
        gridEl.appendChild(cell);
      }
    }
  }

  function legendSwatchHTML(tm) {
    if (tm.single) { return '<i class="bi ' + tm.single + '"></i>'; }
    return tm.symbol || '';
  }

  function typeLegend() {
    const order = ['SMP', 'CTRL', 'NEG', 'POS', 'SKIP', 'EMPTY'];
    let html = '';
    order.forEach(function (key) {
      const tm = TYPE_META[key];
      html += '<span class="plate-map-legend-item">' +
        '<span class="plate-map-legend-swatch" style="background:' + tm.fill + ';border-color:' + tm.border + ';color:' + tm.text + '">' + legendSwatchHTML(tm) + '</span>' +
        tm.label +
        '</span>';
    });
    return html;
  }

  function render() {
    const m = data.metrics[activeMetric];
    const stats = computeStats(m);

    // value (gradient) legend
    if (legendEl) {
      if (m.mode === 'diverging') {
        const lo = (m.min !== undefined) ? m.min : stats.min;
        const hi = (m.max !== undefined) ? m.max : stats.max;
        legendEl.innerHTML =
          '<span class="small">' + fmt(lo) + '</span>' +
          '<div class="plate-map-scale" style="background:linear-gradient(90deg, hsl(120,65%,45%), hsl(0,65%,45%))"></div>' +
          '<span class="small">' + fmt(hi) + '</span>';
      } else {
        legendEl.innerHTML =
          '<span class="small">' + fmt(stats.min) + '</span>' +
          '<div class="plate-map-scale" style="background:linear-gradient(90deg, hsl(0,65%,92%), hsl(0,65%,35%))"></div>' +
          '<span class="small">' + fmt(stats.max) + '</span>';
      }
    }

    const cells = gridEl.querySelectorAll('.plate-map-cell');

    cells.forEach(function (cell, idx) {
      const wellId = cell.dataset.wellId;
      const w = wellById[wellId];
      const t = effectiveType(w);
      const iconHTML = typeIconHTML(w);

      // type class drives the border colour; heatmap fill is applied below
      cell.className = 'plate-map-cell pmc-' + t.toLowerCase();
      cell.style.background = '';
      cell.style.color = '';

      // fill encodes the value (or stays white), border keeps the type colour
      const bg = w ? colorFor(w, stats, m) : null;
      if (bg) {
        cell.style.background = bg;
        cell.style.color = textColorFor(bg);
      } else {
        cell.style.background = '#ffffff';
        cell.style.color = '#212529';
      }

      cell.innerHTML = iconHTML;
      cell.title = w ? w.id : '';
    });
  }

  // metric selector: buttons for inhib, dropdowns for readouts/concentrations
  const DROPDOWN_GROUPS = [
    { key: 'readout', title: 'Readouts' },
    { key: 'conc', title: 'Concentrations' }
  ];
  const dropdownKeys = DROPDOWN_GROUPS.map(function (g) { return g.key; });
  const buttonMetricNames = metricNames.filter(function (n) { return dropdownKeys.indexOf(data.metrics[n].key) === -1; });
  const metricButtons = {};
  const dropdownState = {};

  function setActive(name) {
    activeMetric = name;
    const m = data.metrics[name];
    Object.keys(metricButtons).forEach(function (n) {
      metricButtons[n].classList.toggle('active', n === name);
    });
    DROPDOWN_GROUPS.forEach(function (g) {
      const st = dropdownState[g.key];
      if (!st) { return; }
      st.toggle.classList.toggle('active', m.key === g.key);
      st.items.forEach(function (item) {
        item.classList.toggle('active', item.dataset.metric === name);
      });
    });
    render();
  }

  buttonMetricNames.forEach(function (name) {
    const btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'btn btn-outline-secondary btn-sm';
    btn.textContent = data.metrics[name].label;
    btn.addEventListener('click', function () { setActive(name); });
    metricButtons[name] = btn;
    metricsEl.appendChild(btn);
  });

  DROPDOWN_GROUPS.forEach(function (g) {
    const names = metricNames.filter(function (n) { return data.metrics[n].key === g.key; });
    if (!names.length) { return; }

    const group = document.createElement('div');
    group.className = 'btn-group btn-group-sm';
    const toggle = document.createElement('button');
    toggle.type = 'button';
    toggle.className = 'btn btn-outline-secondary btn-sm dropdown-toggle';
    toggle.setAttribute('data-bs-toggle', 'dropdown');
    toggle.setAttribute('aria-expanded', 'false');
    toggle.textContent = g.title;

    const menu = document.createElement('ul');
    menu.className = 'dropdown-menu';
    const items = [];
    names.forEach(function (name) {
      const li = document.createElement('li');
      const a = document.createElement('a');
      a.className = 'dropdown-item';
      a.href = '#';
      a.dataset.metric = name;
      a.textContent = data.metrics[name].label;
      a.addEventListener('click', function (e) { e.preventDefault(); setActive(name); });
      li.appendChild(a);
      menu.appendChild(li);
      items.push(a);
    });
    group.appendChild(toggle);
    group.appendChild(menu);
    metricsEl.appendChild(group);
    dropdownState[g.key] = { toggle: toggle, items: items };
  });

  // type legend is static and shown across all heatmap tabs (TestPlate only)
  if (typeLegendEl) {
    typeLegendEl.innerHTML = (data.show_type !== false) ? typeLegend() : '';
  }

  // tooltip via event delegation
  gridEl.addEventListener('mousemove', function (e) {
    const cell = e.target.closest ? e.target.closest('.plate-map-cell') : null;
    if (cell) {
      const w = wellById[cell.dataset.wellId];
      if (w) { showTooltip(w, cell); }
      else { hideTooltip(); }
    } else {
      hideTooltip();
    }
  });
  gridEl.addEventListener('mouseleave', hideTooltip);

  // map well data onto cells
  buildGrid();
  (function attachCells() {
    const cells = gridEl.querySelectorAll('.plate-map-cell');
    cells.forEach(function (cell, idx) {
      const r = Math.floor(idx / data.n_cols);
      const c = (idx % data.n_cols) + 1;
      const wellId = data.row_labels[r] + String(c).padStart(2, '0');
      cell.dataset.wellId = wellId;
    });
  })();

  setActive(metricNames[0]);
})();
