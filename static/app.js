document.addEventListener('DOMContentLoaded', () => {
  const filterCatalog = JSON.parse(document.getElementById('filter-catalog').textContent);
  const geneQueryEl = document.getElementById('gene-query');
  const genePlotTypeEl = document.getElementById('gene-plot-type');
  const umapColumnEl = document.getElementById('umap-column');
  const umapPlotEl = document.getElementById('umap-plot');
  const umapImageEl = document.getElementById('umap-image');
  const umapEmptyEl = document.getElementById('umap-empty');
  const umapStatusEl = document.getElementById('umap-status');
  const genePlotEl = document.getElementById('gene-plot');

  const geneDiseaseEl = document.getElementById('gene-disease');
  const geneCellTypeEl = document.getElementById('gene-cell-type');

  const degDiseaseEl = document.getElementById('deg-disease');
  const degCellTypeEl = document.getElementById('deg-cell-type');
  const degRegulationEl = document.getElementById('deg-regulation');
  const degHeadEl = document.getElementById('deg-head');
  const degBodyEl = document.getElementById('deg-body');
  const degStatusEl = document.getElementById('deg-status');
  const volcanoDiseaseEl = document.getElementById('volcano-disease');
  const volcanoImageEl = document.getElementById('volcano-image');
  const volcanoEmptyEl = document.getElementById('volcano-empty');
  const volcanoCaptionEl = document.getElementById('volcano-caption');
  const requestedDisease = new URLSearchParams(window.location.search).get('disease')?.trim() || '';

  const DEG_COLUMNS = [
    'gene',
    'protein',
    'sig',
    'de_pval',
    'de_coef',
    'de_se',
    'dv_pval',
    'dv_coef',
    'dv_se',
    'baseMean',
  ];
  const DEG_ALL_CELL_TYPES = '__all__';
  let degColumns = DEG_COLUMNS.slice();
  const degSort = { column: '', order: '' };

  const populateSelect = (select, values, placeholder) => {
    if (!select) return;
    select.innerHTML = '';
    if (!values || values.length === 0) {
      const empty = document.createElement('option');
      empty.textContent = placeholder;
      empty.value = '';
      select.appendChild(empty);
      return;
    }

    values.forEach((value) => {
      const option = document.createElement('option');
      option.value = value;
      option.textContent = value;
      select.appendChild(option);
    });
  };

  const ALL_FILTER_VALUE = '__all__';
  let lastGenePlot = null;
  let geneRequestId = 0;
  let geneRequestController = null;
  const geneLoadingEl = document.getElementById('gene-loading');
  const geneDownloadEl = document.getElementById('gene-download');

  const prependAllOption = (select, label) => {
    if (!select) return;
    const option = document.createElement('option');
    option.value = ALL_FILTER_VALUE;
    option.textContent = label;
    select.prepend(option);
    option.selected = true;
  };

  const selectRequestedDisease = (select) => {
    if (!requestedDisease || !select) return false;
    const hasRequestedDisease = Array.from(select.options).some((option) => option.value === requestedDisease);
    if (!hasRequestedDisease) return false;
    select.value = requestedDisease;
    return true;
  };

  // An empty value means "no filter" server-side, which is what the All option selects.
  const filterValues = (select) => {
    if (!select || !select.selectedOptions) return '';
    const values = Array.from(select.selectedOptions).map((option) => option.value).filter(Boolean);
    if (values.length === 0 || values.includes(ALL_FILTER_VALUE)) return '';
    return values.join(',');
  };

  const metadataLabel = (column) => column.replaceAll('_', ' ');

  umapImageEl.addEventListener('error', () => {
    umapImageEl.hidden = true;
    umapEmptyEl.hidden = false;
    umapEmptyEl.textContent = 'The pre-generated UMAP image could not be loaded.';
    umapStatusEl.textContent = `Unable to show ${metadataLabel(umapColumnEl.value)}`;
    umapPlotEl.classList.add('empty');
  });
  umapImageEl.addEventListener('load', () => {
    umapImageEl.hidden = false;
    umapEmptyEl.hidden = true;
    umapStatusEl.textContent = `Showing ${metadataLabel(umapColumnEl.value)}`;
    umapPlotEl.classList.remove('empty');
  });

  const updateUmapPlot = () => {
    const column = umapColumnEl?.value;
    if (!column) return;
    const label = metadataLabel(column);
    umapStatusEl.textContent = `Loading ${label}`;
    umapImageEl.alt = `UMAP colored by ${label}`;
    umapImageEl.src = `/umap_image/${encodeURIComponent(column)}.png`;
  };

  const updateGenePlot = async () => {
    const requestId = ++geneRequestId;
    geneRequestController?.abort();
    geneRequestController = null;
    lastGenePlot = null;
    geneDownloadEl.disabled = true;
    geneLoadingEl.hidden = true;
    genePlotEl.setAttribute('aria-busy', 'false');
    const gene = geneQueryEl.value.trim();
    if (!gene) {
      genePlotEl.classList.add('empty');
      genePlotEl.textContent = 'Enter a gene to generate an expression plot.';
      return;
    }

    const disease = filterValues(geneDiseaseEl);
    const cellType = filterValues(geneCellTypeEl);
    const plotType = genePlotTypeEl?.value || 'histogram';
    const controller = new AbortController();
    geneRequestController = controller;
    geneLoadingEl.hidden = false;
    genePlotEl.setAttribute('aria-busy', 'true');
    genePlotEl.classList.add('empty');
    genePlotEl.textContent = 'Please wait while the expression plot is generated.';

    try {
      const params = new URLSearchParams({ gene, disease, cell_type: cellType, plot_type: plotType });
      const response = await fetch(`/api/gene_expression?${params.toString()}`, { signal: controller.signal });
      const payload = await response.json();
      if (requestId !== geneRequestId) return;
      if (!response.ok || payload.error) {
        throw new Error(payload.error || 'Gene expression plot could not be generated.');
      }
      const image = new Image();
      image.alt = `Expression for ${payload.gene}`;
      image.src = payload.image_base64;
      await image.decode();
      if (requestId !== geneRequestId) return;
      genePlotEl.classList.remove('empty');
      genePlotEl.replaceChildren(image);
      lastGenePlot = { gene: payload.gene, plotType, image: payload.image_base64 };
      geneDownloadEl.disabled = false;
    } catch (error) {
      if (requestId !== geneRequestId) return;
      genePlotEl.classList.add('empty');
      genePlotEl.textContent = error instanceof Error ? error.message : 'Gene expression plot could not be generated.';
      lastGenePlot = null;
    } finally {
      if (requestId === geneRequestId) {
        geneRequestController = null;
        geneLoadingEl.hidden = true;
        genePlotEl.setAttribute('aria-busy', 'false');
      }
    }
  };

  const downloadGenePlot = () => {
    if (!lastGenePlot) return;
    const link = document.createElement('a');
    link.href = lastGenePlot.image;
    link.download = `${lastGenePlot.gene}_${lastGenePlot.plotType}_expression.png`;
    link.click();
  };

  const selectedValues = (select) =>
    select && select.selectedOptions ? Array.from(select.selectedOptions).map((option) => option.value).filter(Boolean).join(',') : '';

  const renderDegHead = () => {
    degHeadEl.replaceChildren();
    degColumns.forEach((column) => {
      const th = document.createElement('th');
      th.className = 'sortable';
      th.dataset.column = column;
      th.tabIndex = 0;
      const arrow = degSort.column === column && degSort.order ? (degSort.order === 'asc' ? ' \u25B2' : ' \u25BC') : ' \u21C5';
      th.textContent = column + arrow;
      th.addEventListener('click', () => cycleDegSort(column));
      th.addEventListener('keydown', (event) => {
        if (event.key === 'Enter' || event.key === ' ') {
          event.preventDefault();
          cycleDegSort(column);
        }
      });
      degHeadEl.appendChild(th);
    });
  };

  const formatDegValue = (column, value) => {
    if (value === null || value === undefined || value === '') return '';
    if (typeof value === 'number' && Number.isFinite(value)) {
      const magnitude = Math.abs(value);
      if (magnitude !== 0 && (magnitude < 1e-3 || magnitude >= 1e5)) return value.toExponential(3);
      return String(Number(value.toFixed(4)));
    }
    const text = String(value);
    return column === 'protein' && text.length > 120 ? `${text.slice(0, 120)}...` : text;
  };

  const renderDegRows = (rows) => {
    degColumns = rows && rows.length > 0 && 'cell_type' in rows[0]
      ? ['gene', 'cell_type', ...DEG_COLUMNS.filter((column) => column !== 'gene')]
      : DEG_COLUMNS.slice();
    renderDegHead();
    degBodyEl.replaceChildren();
    if (!rows || rows.length === 0) {
      const row = document.createElement('tr');
      const cell = document.createElement('td');
      cell.colSpan = degColumns.length;
      cell.textContent = 'No differentially expressed genes for these filters.';
      row.appendChild(cell);
      degBodyEl.appendChild(row);
      return;
    }

    rows.forEach((record) => {
      const row = document.createElement('tr');
      degColumns.forEach((column) => {
        const cell = document.createElement('td');
        cell.textContent = formatDegValue(column, record[column]);
        if (column === 'sig') cell.className = `sig-${String(record[column] ?? '').toLowerCase()}`;
        row.appendChild(cell);
      });
      degBodyEl.appendChild(row);
    });
  };

  const degQueryParams = () =>
    new URLSearchParams({
      disease: selectedValues(degDiseaseEl),
      cell_type: selectedValues(degCellTypeEl),
      regulation: degRegulationEl?.value || '',
      sort_by: degSort.order ? degSort.column : '',
      sort_order: degSort.order || 'desc',
    });

  const loadDegTable = async () => {
    if (!degHeadEl || !degBodyEl) return;

    const disease = selectedValues(degDiseaseEl);
    if (!disease) {
      renderDegRows([]);
      degStatusEl.textContent = 'Select at least one disease to load DEG results.';
      return;
    }

    degStatusEl.textContent = 'Loading DEG results...';

    const params = degQueryParams();

    try {
      const response = await fetch(`/api/deg_table?${params.toString()}`);
      const payload = await response.json();
      renderDegRows(payload.rows || []);
      degStatusEl.textContent = `${payload.count ?? 0} rows shown`;
    } catch (error) {
      renderDegRows([]);
      degStatusEl.textContent = 'DEG results could not be loaded.';
    }
  };

  const cycleDegSort = (column) => {
    if (degSort.column !== column) {
      degSort.column = column;
      degSort.order = 'asc';
    } else if (degSort.order === 'asc') {
      degSort.order = 'desc';
    } else if (degSort.order === 'desc') {
      degSort.column = '';
      degSort.order = '';
    } else {
      degSort.order = 'asc';
    }
    loadDegTable();
  };

  const refreshDegOptions = () => {
    if (!degDiseaseEl || !degCellTypeEl) return;
    populateSelect(degDiseaseEl, filterCatalog.deg.diseases, 'No diseases');
    populateSelect(degCellTypeEl, filterCatalog.deg.cell_types, 'No cell types');
    const diseasePlaceholder = document.createElement('option');
    diseasePlaceholder.value = '';
    diseasePlaceholder.textContent = 'Select a disease';
    degDiseaseEl.prepend(diseasePlaceholder);
    const allOption = document.createElement('option');
    allOption.value = DEG_ALL_CELL_TYPES;
    allOption.textContent = 'All cell types';
    degCellTypeEl.prepend(allOption);
    degDiseaseEl.value = '';
    degCellTypeEl.value = DEG_ALL_CELL_TYPES;
    selectRequestedDisease(degDiseaseEl);
  };

  const loadVolcanoCatalog = () => {
    if (!volcanoDiseaseEl) return;
    populateSelect(volcanoDiseaseEl, filterCatalog.volcano_diseases, 'No volcano plots available');
    const diseases = filterCatalog.volcano_diseases;
    if (diseases.length === 0) {
      volcanoCaptionEl.textContent = 'No volcano plots available';
      volcanoEmptyEl.textContent = 'Volcano plots are unavailable.';
      return;
    }
    volcanoDiseaseEl.value = selectRequestedDisease(volcanoDiseaseEl) ? requestedDisease : diseases[0];
    renderVolcanoPlot();
  };

  const renderVolcanoPlot = () => {
    const disease = volcanoDiseaseEl?.value || '';
    if (!disease) return;
    volcanoImageEl.src = `/volcano_image/${encodeURIComponent(disease)}.png`;
    volcanoImageEl.alt = `Volcano plot for ${disease}`;
    volcanoImageEl.hidden = false;
    volcanoEmptyEl.hidden = true;
    volcanoCaptionEl.textContent = disease;
  };

  const bind = () => {
    document.getElementById('gene-plot-button')?.addEventListener('click', updateGenePlot);
    genePlotTypeEl?.addEventListener('change', updateGenePlot);
    geneDiseaseEl?.addEventListener('change', updateGenePlot);
    geneCellTypeEl?.addEventListener('change', updateGenePlot);
    document.getElementById('gene-download')?.addEventListener('click', downloadGenePlot);
    document.getElementById('gene-reset')?.addEventListener('click', () => {
      if (geneDiseaseEl) geneDiseaseEl.value = ALL_FILTER_VALUE;
      if (geneCellTypeEl) geneCellTypeEl.value = ALL_FILTER_VALUE;
      geneQueryEl.value = '';
      updateGenePlot();
    });
    geneQueryEl.addEventListener('keydown', (event) => {
      if (event.key === 'Enter') updateGenePlot();
    });

    document.getElementById('deg-apply')?.addEventListener('click', loadDegTable);
    document.getElementById('deg-download')?.addEventListener('click', () => {
      if (!selectedValues(degDiseaseEl)) {
        degStatusEl.textContent = 'Select at least one disease before downloading.';
        return;
      }
      window.location.href = `/api/deg_download?${degQueryParams().toString()}`;
    });
    document.getElementById('deg-reset')?.addEventListener('click', () => {
      if (degDiseaseEl) degDiseaseEl.value = '';
      if (degCellTypeEl) degCellTypeEl.value = DEG_ALL_CELL_TYPES;
      if (degRegulationEl) degRegulationEl.value = '';
      degSort.column = '';
      degSort.order = '';
      loadDegTable();
    });
    degRegulationEl?.addEventListener('change', loadDegTable);
    degDiseaseEl?.addEventListener('change', loadDegTable);
    degCellTypeEl?.addEventListener('change', loadDegTable);
    volcanoDiseaseEl?.addEventListener('change', renderVolcanoPlot);
    umapColumnEl?.addEventListener('change', updateUmapPlot);
  };

  populateSelect(geneCellTypeEl, filterCatalog.metadata.cell_types, 'No cell types');
  prependAllOption(geneCellTypeEl, 'All cell types');
  populateSelect(geneDiseaseEl, filterCatalog.metadata.diseases, 'No diseases');
  prependAllOption(geneDiseaseEl, 'All diseases');
  if (selectRequestedDisease(geneDiseaseEl)) updateGenePlot();
  refreshDegOptions();
  if (requestedDisease && degDiseaseEl?.value === requestedDisease) {
    loadDegTable();
  } else {
    renderDegRows([]);
    degStatusEl.textContent = 'Select a disease to load DEGs.';
  }
  loadVolcanoCatalog();
  bind();
  updateGenePlot();
});
