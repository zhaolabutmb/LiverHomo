document.addEventListener('DOMContentLoaded', () => {
  const filterCatalog = JSON.parse(document.getElementById('filter-catalog').textContent);
  const graphCanvas = document.getElementById('graph-canvas');
  const diseaseSelect = document.getElementById('graph-disease');
  const cellTypeSelect = document.getElementById('graph-cell-type');
  const geneListEl = document.getElementById('gene-list');
  const graphLayoutSelect = document.getElementById('graph-layout');
  const graphCaption = document.getElementById('graph-caption');
  const graphDownloadSvg = document.getElementById('graph-download-svg');
  const graphDownloadSource = document.getElementById('graph-download-source');
  const graphZoomIn = document.getElementById('graph-zoom-in');
  const graphZoomOut = document.getElementById('graph-zoom-out');
  const graphZoomReset = document.getElementById('graph-zoom-reset');
  const egoDiseaseSelect = document.getElementById('ego-disease');
  const egoCellTypeSelect = document.getElementById('ego-cell-type');
  const egoImage = document.getElementById('ego-image');
  const egoEmpty = document.getElementById('ego-empty');
  const egoCaption = document.getElementById('ego-caption');
  const interactionModeSelect = document.getElementById('interaction-mode');
  const interactionSelect = document.getElementById('interaction-disease');
  const interactionDiseaseLabel = document.getElementById('interaction-disease-label');
  const interactionCellTypeASelect = document.getElementById('interaction-cell-type-a');
  const interactionCellTypeALabel = document.getElementById('interaction-cell-type-a-label');
  const interactionCellTypeBSelect = document.getElementById('interaction-cell-type-b');
  const interactionCellTypeBLabel = document.getElementById('interaction-cell-type-b-label');
  const interactionImage = document.getElementById('interaction-image');
  const interactionEmpty = document.getElementById('interaction-empty');
  const interactionCaption = document.getElementById('interaction-caption');
  let interactionCellTypes = [];
  let currentGraph = null;
  let graphViewport = null;
  let graphView = { x: 0, y: 0, scale: 1 };
  let panStart = null;
  let graphRequestId = 0;
  let graphLoadTimer = null;
  let graphRequestController = null;
  const graphLoadingEl = document.getElementById('graph-loading');

  const applyGraphView = () => {
    if (!graphViewport) return;
    graphViewport.setAttribute('transform', `translate(${graphView.x} ${graphView.y}) scale(${graphView.scale})`);
  };

  const resetGraphView = () => {
    graphView = { x: 0, y: 0, scale: 1 };
    applyGraphView();
  };

  const zoomGraph = (factor, center = { x: 440, y: 320 }) => {
    const nextScale = Math.max(0.35, Math.min(6, graphView.scale * factor));
    const ratio = nextScale / graphView.scale;
    graphView = {
      x: center.x - (center.x - graphView.x) * ratio,
      y: center.y - (center.y - graphView.y) * ratio,
      scale: nextScale,
    };
    applyGraphView();
  };

  const svgPointFromEvent = (event) => {
    const rect = graphCanvas.getBoundingClientRect();
    return {
      x: ((event.clientX - rect.left) / rect.width) * 880,
      y: ((event.clientY - rect.top) / rect.height) * 640,
    };
  };

  const startPan = (event, pointerId = 'mouse') => {
    if (!graphViewport) return;
    panStart = { pointerId, x: event.clientX, y: event.clientY, viewX: graphView.x, viewY: graphView.y };
  };

  const movePan = (event, pointerId = 'mouse') => {
    if (!panStart || panStart.pointerId !== pointerId) return;
    const rect = graphCanvas.getBoundingClientRect();
    graphView.x = panStart.viewX + ((event.clientX - panStart.x) / rect.width) * 880;
    graphView.y = panStart.viewY + ((event.clientY - panStart.y) / rect.height) * 640;
    applyGraphView();
  };

  const endPan = (pointerId = 'mouse') => {
    if (panStart && panStart.pointerId === pointerId) panStart = null;
  };

  const setGraphCaption = (text) => {
    if (graphCaption) graphCaption.textContent = text;
  };

  const setDownloadState = () => {
    const hasGraph = currentGraph && currentGraph.nodes && currentGraph.nodes.length > 0;
    if (graphDownloadSvg) graphDownloadSvg.disabled = !hasGraph;
    if (graphDownloadSource) graphDownloadSource.disabled = !hasGraph || !currentGraph.graph || !currentGraph.graph.file;
  };

  const loadGraphSummary = () => {
    const { diseases, cell_types: cellTypes } = filterCatalog.graph;
    populateSelect(diseaseSelect, diseases, 'Select disease');
    populateSelect(cellTypeSelect, cellTypes, 'Select cell type');
    if (diseases.length > 0) {
      diseaseSelect.value = diseases[0];
    }
    if (cellTypes.length > 0) {
      cellTypeSelect.value = cellTypes.includes('CD8_EM') ? 'CD8_EM' : cellTypes[0];
    }
  };

  const populateSelect = (select, values, emptyLabel = '') => {
    if (!select) return;
    select.innerHTML = '';
    if (emptyLabel) {
      const emptyOption = document.createElement('option');
      emptyOption.value = '';
      emptyOption.textContent = emptyLabel;
      select.appendChild(emptyOption);
    }
    if (!values || values.length === 0) {
      const option = document.createElement('option');
      option.value = '';
      option.textContent = 'No entries';
      select.appendChild(option);
      return;
    }
    values.forEach((value) => {
      const option = document.createElement('option');
      option.value = value;
      option.textContent = value;
      select.appendChild(option);
    });
  };

  const renderGraph = (graph) => {
    if (!graphCanvas) return;
    graphCanvas.innerHTML = '';
    graphViewport = null;
    panStart = null;
    currentGraph = graph;
    const loading = Boolean(graph && graph.loading);
    graphLoadingEl.hidden = !loading;
    graphCanvas.setAttribute('aria-busy', String(loading));

    if (graph && graph.error) {
      graphCanvas.innerHTML = `<text x="50%" y="50%" text-anchor="middle" fill="#9a3412" font-size="16">${graph.error}</text>`;
      setGraphCaption('Graph unavailable');
      setDownloadState();
      return;
    }

    if (graph && graph.loading) {
      graphCanvas.innerHTML = '<text x="50%" y="50%" text-anchor="middle" fill="#5a7287" font-size="18">Loading graph data...</text>';
      setGraphCaption('Loading graph data');
      setDownloadState();
      return;
    }

    if (!graph || !graph.nodes || graph.nodes.length === 0) {
      graphCanvas.innerHTML = '<text x="50%" y="50%" text-anchor="middle" fill="#5a7287" font-size="18">Gene not found in current graph filters.</text>';
      setGraphCaption('No matching graph data');
      setDownloadState();
      return;
    }

    const width = 880;
    const height = 640;
    const marginX = 64;
    const marginY = 46;
    const plotWidth = width - marginX * 2;
    const plotHeight = height - marginY * 2;
    const nodes = graph.nodes;
    const edges = graph.edges || [];

    const positions = {};
    nodes.forEach((node) => {
      positions[node.id] = { x: node.x, y: node.y };
    });

    const svgNS = 'http://www.w3.org/2000/svg';
    graphViewport = document.createElementNS(svgNS, 'g');
    graphViewport.setAttribute('class', 'graph-viewport');
    const edgeLayer = document.createElementNS(svgNS, 'g');
    const nodeLayer = document.createElementNS(svgNS, 'g');
    const labelLayer = document.createElementNS(svgNS, 'g');

    edges.forEach((edge) => {
      const a = positions[edge.source];
      const b = positions[edge.target];
      if (!a || !b) return;
      const line = document.createElementNS(svgNS, 'line');
      line.setAttribute('x1', a.x * plotWidth + marginX);
      line.setAttribute('y1', (1 - a.y) * plotHeight + marginY);
      line.setAttribute('x2', b.x * plotWidth + marginX);
      line.setAttribute('y2', (1 - b.y) * plotHeight + marginY);
      line.setAttribute('class', 'edge');
      edgeLayer.appendChild(line);
    });

    nodes.forEach((node) => {
      const cx = node.x * plotWidth + marginX;
      const cy = (1 - node.y) * plotHeight + marginY;
      const radius = Math.max(6, Math.min(18, 6 + node.degree * 1.5));
      const circle = document.createElementNS(svgNS, 'circle');
      circle.setAttribute('cx', cx);
      circle.setAttribute('cy', cy);
      circle.setAttribute('r', radius);
      circle.setAttribute('class', 'node');
      circle.setAttribute('tabindex', '0');
      circle.setAttribute('aria-label', `${node.id}, degree ${node.degree}`);
      const title = document.createElementNS(svgNS, 'title');
      title.textContent = `${node.id} degree ${node.degree}`;
      circle.appendChild(title);
      nodeLayer.appendChild(circle);

      const label = document.createElementNS(svgNS, 'text');
      const directionX = cx - width / 2;
      const directionY = cy - height / 2;
      const directionLength = Math.hypot(directionX, directionY) || 1;
      const labelOffset = radius + 10;
      label.setAttribute('x', cx + (directionX / directionLength) * labelOffset);
      label.setAttribute('y', cy + (directionY / directionLength) * labelOffset);
      label.setAttribute('class', 'node-label');
      label.setAttribute('text-anchor', directionX >= 0 ? 'start' : 'end');
      label.setAttribute('dominant-baseline', 'middle');
      label.textContent = node.id;
      labelLayer.appendChild(label);
    });

    graphViewport.appendChild(edgeLayer);
    graphViewport.appendChild(nodeLayer);
    graphViewport.appendChild(labelLayer);
    graphCanvas.appendChild(graphViewport);
    resetGraphView();
    const label = graph.graph ? `${graph.graph.disease} / ${graph.graph.cell_type}` : 'Current graph';
    const layoutLabel = graphLayoutSelect.selectedOptions[0]?.textContent || 'Random';
    const total = graph.total_edges || edges.length;
    const shown = total > edges.length
      ? `top ${edges.length.toLocaleString()} of ${total.toLocaleString()} interactions by PrePPI score`
      : `${edges.length.toLocaleString()} interactions`;
    setGraphCaption(`${label} · ${layoutLabel}: ${nodes.length} nodes, ${shown}`);
    setDownloadState();
  };

  const downloadCurrentSvg = () => {
    if (!graphCanvas || !currentGraph || !currentGraph.nodes || currentGraph.nodes.length === 0) return;
    const clone = graphCanvas.cloneNode(true);
    clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    const style = document.createElementNS('http://www.w3.org/2000/svg', 'style');
    style.textContent = `
      .node { fill: #2c7be5; stroke: rgba(9, 52, 95, 0.2); stroke-width: 1.2; }
      .edge { stroke: rgba(68, 116, 160, 0.5); stroke-width: 1.5; }
      .node-label { font-size: 12px; fill: #173a52; font-weight: 700; font-family: sans-serif; }
    `;
    clone.insertBefore(style, clone.firstChild);
    const source = new XMLSerializer().serializeToString(clone);
    const blob = new Blob([source], { type: 'image/svg+xml;charset=utf-8' });
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    const name = currentGraph.graph ? `${currentGraph.graph.disease}_${currentGraph.graph.cell_type}_network.svg` : 'graph_network.svg';
    link.href = url;
    link.download = name.replaceAll(' ', '_');
    document.body.appendChild(link);
    link.click();
    link.remove();
    URL.revokeObjectURL(url);
  };

  const downloadSourceGraph = () => {
    if (!currentGraph || !currentGraph.graph || !currentGraph.graph.file) return;
    const params = new URLSearchParams({ file: currentGraph.graph.file });
    window.location.href = `/api/graph_download?${params.toString()}`;
  };

  const loadGraph = async () => {
    const requestId = ++graphRequestId;
    graphRequestController?.abort();
    const controller = new AbortController();
    graphRequestController = controller;
    window.clearTimeout(graphLoadTimer);
    const disease = diseaseSelect.value;
    const cellType = Array.from(cellTypeSelect.selectedOptions).map((option) => option.value).join(',');
    const genes = geneListEl.value.trim();
    const layout = graphLayoutSelect.value;
    const params = new URLSearchParams({ disease, cell_type: cellType, gene_list: genes, layout, limit: '50' });
    renderGraph({ loading: true });

    const fetchGraph = async () => {
      if (requestId !== graphRequestId) return;
      try {
        const response = await fetch(`/api/graph_network?${params.toString()}`, { signal: controller.signal });
        const payload = await response.json();
        if (requestId !== graphRequestId) return;
        if (!response.ok || payload.error) {
          renderGraph({ error: payload.error || 'Graph data could not be loaded.' });
          return;
        }
        if (payload.loading) {
          renderGraph(payload);
          graphLoadTimer = window.setTimeout(fetchGraph, 750);
          return;
        }
        renderGraph(payload);
      } catch (error) {
        if (requestId !== graphRequestId) return;
        renderGraph({ error: 'Graph data could not be loaded. Please try again.' });
      }
    };
    await fetchGraph();
  };

  const renderEgo = () => {
    const disease = egoDiseaseSelect.value;
    const cellType = egoCellTypeSelect.value;
    if (!disease || !cellType) {
      egoImage.removeAttribute('src');
      egoImage.hidden = true;
      egoEmpty.hidden = false;
      egoCaption.textContent = 'Select a disease and cell type';
      return;
    }
    egoImage.src = `/ego_image/${encodeURIComponent(disease)}/${encodeURIComponent(cellType)}.png`;
    egoImage.hidden = false;
    egoEmpty.hidden = true;
    egoCaption.textContent = `${disease} / ${cellType.replaceAll('_', ' ')}`;
  };

  const updateInteractionControls = () => {
    const mode = interactionModeSelect.value;
    const isHistogram = mode === 'histogram';
    interactionDiseaseLabel.hidden = isHistogram;
    interactionCellTypeALabel.hidden = !isHistogram;
    interactionCellTypeBLabel.hidden = !isHistogram;
  };

  const renderInteraction = async () => {
    const mode = interactionModeSelect.value;

    if (mode === 'histogram') {
      const cellTypeA = interactionCellTypeASelect.value;
      const cellTypeB = interactionCellTypeBSelect.value;
      if (!cellTypeA || !cellTypeB) {
        interactionImage.removeAttribute('src');
        interactionImage.hidden = true;
        interactionEmpty.hidden = false;
        interactionCaption.textContent = 'Select two cell types';
        return;
      }
      try {
        const response = await fetch(
          `/api/cell_interaction_pair?cell_type_a=${encodeURIComponent(cellTypeA)}&cell_type_b=${encodeURIComponent(cellTypeB)}`
        );
        const payload = await response.json();
        if (!response.ok) throw new Error(payload.error || 'Interaction data unavailable');
        interactionImage.src = payload.image_base64;
        interactionImage.hidden = false;
        interactionEmpty.hidden = true;
        interactionCaption.textContent = `${cellTypeA} \u2194 ${cellTypeB} by disease`;
      } catch (error) {
        interactionImage.removeAttribute('src');
        interactionImage.hidden = true;
        interactionEmpty.hidden = false;
        interactionCaption.textContent = 'Interaction data unavailable';
      }
      return;
    }

    const disease = interactionSelect.value;
    if (!disease) {
      interactionImage.removeAttribute('src');
      interactionImage.hidden = true;
      interactionEmpty.hidden = false;
      interactionCaption.textContent = 'Select a disease';
      return;
    }
    interactionImage.src = `/cell_interaction/${encodeURIComponent(disease)}.png?kind=${encodeURIComponent(mode)}`;
    interactionImage.hidden = false;
    interactionEmpty.hidden = true;
    interactionCaption.textContent = `${disease} (${mode === 'chord' ? 'chord plot' : 'heatmap'})`;
  };

  const downloadInteraction = () => {
    const mode = interactionModeSelect.value;
    if (mode === 'histogram') {
      if (!interactionImage.src) return;
      const link = document.createElement('a');
      link.href = interactionImage.src;
      link.download = `${interactionCellTypeASelect.value}_${interactionCellTypeBSelect.value}_interactions.png`;
      link.click();
      return;
    }
    const disease = interactionSelect.value;
    if (!disease) return;
    window.location.href = `/cell_interaction/${encodeURIComponent(disease)}.png?kind=${encodeURIComponent(mode)}&download=1`;
  };

  const loadInteractionCatalog = () => {
    if (!interactionSelect) return;
    const { diseases, cell_types: cellTypes } = filterCatalog.cell_interactions;
    interactionCellTypes = cellTypes;
    populateSelect(interactionSelect, diseases);
    populateSelect(interactionCellTypeASelect, interactionCellTypes);
    populateSelect(interactionCellTypeBSelect, interactionCellTypes);
    interactionModeSelect.addEventListener('change', () => {
      updateInteractionControls();
      renderInteraction();
    });
    interactionSelect.addEventListener('change', renderInteraction);
    interactionCellTypeASelect.addEventListener('change', renderInteraction);
    interactionCellTypeBSelect.addEventListener('change', renderInteraction);
    document.getElementById('interaction-view-btn').addEventListener('click', renderInteraction);
    document.getElementById('interaction-download').addEventListener('click', downloadInteraction);
    updateInteractionControls();
    if (diseases.length > 0) {
      interactionSelect.value = diseases[0];
      renderInteraction();
    }
  };


  const loadEgoCatalog = () => {
    const graphs = filterCatalog.ego;
    const diseases = [...new Set(graphs.map((item) => item.disease))];
    populateSelect(egoDiseaseSelect, diseases);
    egoDiseaseSelect.addEventListener('change', () => {
      const cellTypes = graphs
        .filter((item) => item.disease === egoDiseaseSelect.value)
        .map((item) => item.cell_type);
      populateSelect(egoCellTypeSelect, cellTypes);
      renderEgo();
    });
    egoCellTypeSelect.addEventListener('change', renderEgo);
    document.getElementById('ego-view').addEventListener('click', renderEgo);
    if (diseases.length > 0) {
      egoDiseaseSelect.value = diseases[0];
      egoDiseaseSelect.dispatchEvent(new Event('change'));
    }
  };

  document.getElementById('graph-apply').addEventListener('click', loadGraph);
  graphLayoutSelect.addEventListener('change', loadGraph);
  if (graphDownloadSvg) graphDownloadSvg.addEventListener('click', downloadCurrentSvg);
  if (graphDownloadSource) graphDownloadSource.addEventListener('click', downloadSourceGraph);
  if (graphZoomIn) graphZoomIn.addEventListener('click', () => zoomGraph(1.25));
  if (graphZoomOut) graphZoomOut.addEventListener('click', () => zoomGraph(0.8));
  if (graphZoomReset) graphZoomReset.addEventListener('click', resetGraphView);
  if (graphCanvas) {
    graphCanvas.addEventListener('wheel', (event) => {
      event.preventDefault();
      zoomGraph(event.deltaY < 0 ? 1.12 : 0.88, svgPointFromEvent(event));
    });
    graphCanvas.addEventListener('pointerdown', (event) => {
      startPan(event, event.pointerId);
    });
    graphCanvas.addEventListener('pointermove', (event) => {
      movePan(event, event.pointerId);
    });
    graphCanvas.addEventListener('pointerup', (event) => {
      endPan(event.pointerId);
    });
    graphCanvas.addEventListener('pointercancel', () => {
      panStart = null;
    });
    graphCanvas.addEventListener('mousedown', (event) => {
      startPan(event);
    });
    graphCanvas.addEventListener('mousemove', (event) => {
      if (event.buttons !== 1) return;
      movePan(event);
    });
    document.addEventListener('pointermove', (event) => {
      movePan(event, event.pointerId);
    });
    document.addEventListener('pointerup', (event) => {
      endPan(event.pointerId);
    });
    document.addEventListener('mouseup', () => {
      endPan();
    });
  }
  document.getElementById('graph-reset').addEventListener('click', () => {
    if (diseaseSelect) diseaseSelect.value = '';
    if (cellTypeSelect) cellTypeSelect.selectedIndex = -1;
    geneListEl.value = '';
    graphLayoutSelect.value = 'random';
    loadGraph();
  });

  loadGraphSummary();
  loadGraph();
  loadEgoCatalog();
  loadInteractionCatalog();
});
