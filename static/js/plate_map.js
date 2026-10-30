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
  const statusEl = document.getElementById('plate-map-status');
  const highlightEl = document.getElementById('plate-map-highlight');
  const tooltipEl = document.getElementById('plate-map-tooltip');

  if (!gridEl) { return; }

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
  const showValue = (data.n_wells || wells.length) <= 96;

  const wellById = {};
  wells.forEach(function (w) { wellById[w.id] = w; });

  // ---------- helpers ----------
  function metricValue(w, metric) {
    if (metric === 'inhibition') {
      const v = w.inhibition;
      return (v === null || v === undefined || !isFinite(Number(v))) ? null : Number(v);
    }
    if (metric === 'conc') {
      const arr = w.conc || [];
      return arr.length ? Number(arr[0]) : null;
    }
    if (metric === 'readout') {
      const arr = w.readout || [];
      return arr.length ? Number(arr[0]) : null;
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

  function computeStats(metric) {
    const vals = [];
    wells.forEach(function (w) {
      const v = metricValue(w, metric);
      if (v !== null) { vals.push(v); }
    });
    if (vals.length === 0) { return { min: 0, max: 1 }; }
    return { min: Math.min.apply(null, vals), max: Math.max.apply(null, vals) };
  }

  function colorFor(w, metric, stats, m) {
    const v = metricValue(w, metric);
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
    const m = bg.match(/hsl\(0,\s*65%,\s*([\d.]+)%\)/);
    if (m) { return Number(m[1]) > 65 ? '#212529' : '#ffffff'; }
    return '#212529';
  }

  function wellFlags(w) {
    const flags = [];
    if (w.is_negcontrol) { flags.push('NEG'); }
    if (w.is_poscontrol) { flags.push('POS'); }
    if (w.is_control && !w.is_poscontrol && !w.is_negcontrol) { flags.push('CTRL'); }
    if (w.is_sample) { flags.push('SMP'); }
    if (w.is_skip) { flags.push('SKIP'); }
    return flags;
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
      parts.push('<div class="plate-map-tt-row"><span>Compound Batch</span><b>' + items.join('<br>') + '</b></div>');
    } else {
      parts.push('<div class="plate-map-tt-row"><span>Compound Batch</span><b>—</b></div>');
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

  function render() {
    const m = data.metrics[activeMetric];
    const stats = computeStats(activeMetric);

    // legend
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
    const highlight = highlightEl ? highlightEl.checked : true;

    cells.forEach(function (cell, idx) {
      const wellId = cell.dataset.wellId;
      const w = wellById[wellId];
      const bg = w ? colorFor(w, activeMetric, stats, m) : null;

      cell.classList.toggle('has-cmp', !!(w && w.cmpbatch && w.cmpbatch.length));
      cell.classList.toggle('dim', !!(highlight && w && (!w.cmpbatch || !w.cmpbatch.length)));

      if (!bg) {
        cell.style.background = 'transparent';
        cell.style.color = '#212529';
        cell.textContent = '';
        cell.title = w ? w.id : '';
        return;
      }

      cell.style.background = bg;
      cell.style.color = textColorFor(bg);
      cell.textContent = showValue ? fmt(metricValue(w, activeMetric)) : '';
      cell.title = '';
    });
  }

  // metric selector buttons
  metricNames.forEach(function (name) {
    const btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'btn btn-outline-secondary btn-sm';
    btn.textContent = data.metrics[name].label;
    btn.addEventListener('click', function () {
      activeMetric = name;
      metricsEl.querySelectorAll('button').forEach(function (b) {
        b.classList.toggle('active', b === btn);
      });
      render();
    });
    if (name === activeMetric) { btn.classList.add('active'); }
    metricsEl.appendChild(btn);
  });

  // bind highlight toggle
  if (highlightEl) {
    highlightEl.addEventListener('change', render);
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

  render();
})();
