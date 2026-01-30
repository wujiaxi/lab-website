/**
 * H5AD Spatial Data Explorer
 * Interactive visualization of spatial transcriptomics data
 */

// Configuration
const DATA_BASE_URL = '../data';

// State
const state = {
    datasets: [],
    currentDataset: null,
    metadata: null,
    expression: null,
    selectedGene: null,
    colorBy: 'cluster',
    markerSize: 3,
    opacity: 0.8,
    highlightedCategories: new Set(),
};

// Color palettes
const CATEGORICAL_COLORS = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
    '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5',
    '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5',
    '#393b79', '#5254a3', '#6b6ecf', '#9c9ede', '#637939',
    '#8ca252', '#b5cf6b', '#cedb9c', '#8c6d31', '#bd9e39',
    '#e7ba52', '#e7cb94', '#843c39', '#ad494a', '#d6616b',
    '#e7969c', '#7b4173', '#a55194', '#ce6dbd', '#de9ed6'
];

const EXPRESSION_COLORSCALE = [
    [0, '#e0e0e0'],
    [0.25, '#ffeda0'],
    [0.5, '#feb24c'],
    [0.75, '#f03b20'],
    [1, '#bd0026']
];

// DOM Elements
const elements = {
    datasetSelect: document.getElementById('dataset-select'),
    colorBySelect: document.getElementById('color-by-select'),
    geneSearch: document.getElementById('gene-search'),
    geneAutocomplete: document.getElementById('gene-autocomplete'),
    markerSize: document.getElementById('marker-size'),
    markerSizeValue: document.getElementById('marker-size-value'),
    opacity: document.getElementById('opacity'),
    opacityValue: document.getElementById('opacity-value'),
    plot: document.getElementById('plot'),
    loading: document.getElementById('loading'),
    nCells: document.getElementById('n-cells'),
    nGenes: document.getElementById('n-genes'),
    legendTitle: document.getElementById('legend-title'),
    legendList: document.getElementById('legend-list'),
    btnResetZoom: document.getElementById('btn-reset-zoom'),
    btnDownload: document.getElementById('btn-download'),
};

// Initialize
async function init() {
    showLoading(true);

    try {
        // Load dataset index
        const response = await fetch(`${DATA_BASE_URL}/datasets.json`);
        state.datasets = (await response.json()).datasets;

        // Populate dataset dropdown
        populateDatasetDropdown();

        // Load first dataset
        if (state.datasets.length > 0) {
            await loadDataset(state.datasets[0].name);
        }

        // Set up event listeners
        setupEventListeners();

    } catch (error) {
        console.error('Failed to initialize:', error);
        showError('Failed to load data. Please check the console for details.');
    }
}

function populateDatasetDropdown() {
    elements.datasetSelect.innerHTML = state.datasets.map(ds =>
        `<option value="${ds.name}">${ds.display_name || ds.name}</option>`
    ).join('');
}

async function loadDataset(datasetName) {
    showLoading(true, 'Loading metadata...');

    try {
        // Load metadata
        const metaResponse = await fetch(`${DATA_BASE_URL}/${datasetName}_metadata.json`);
        state.metadata = await metaResponse.json();
        state.currentDataset = datasetName;

        // Update stats
        elements.nCells.textContent = state.metadata.n_cells.toLocaleString();
        elements.nGenes.textContent = state.metadata.n_genes.toLocaleString();

        // Update color by options based on available metadata
        updateColorByOptions();

        // Load expression data
        showLoading(true, 'Loading gene expression...');
        const exprResponse = await fetch(`${DATA_BASE_URL}/${datasetName}_expression.json`);
        state.expression = await exprResponse.json();

        // Clear gene selection
        state.selectedGene = null;
        elements.geneSearch.value = '';

        // Render plot
        renderPlot();

    } catch (error) {
        console.error('Failed to load dataset:', error);
        showError(`Failed to load dataset: ${datasetName}`);
    }

    showLoading(false);
}

function updateColorByOptions() {
    const cols = state.metadata.metadata_columns || [];
    elements.colorBySelect.innerHTML = cols.map(col =>
        `<option value="${col}">${formatColumnName(col)}</option>`
    ).join('');

    // Default to cluster if available
    if (cols.includes('cluster')) {
        elements.colorBySelect.value = 'cluster';
        state.colorBy = 'cluster';
    } else if (cols.length > 0) {
        state.colorBy = cols[0];
    }
}

function formatColumnName(col) {
    return col.split('.').map(s =>
        s.charAt(0).toUpperCase() + s.slice(1)
    ).join(' ');
}

function renderPlot() {
    const { coordinates, metadata } = state.metadata;
    const colorByData = metadata[state.colorBy] || [];

    let data, layout;

    if (state.selectedGene && state.expression?.genes?.[state.selectedGene]) {
        // Gene expression mode
        const expr = state.expression.genes[state.selectedGene];

        data = [{
            x: coordinates.x,
            y: coordinates.y,
            mode: 'markers',
            type: 'scattergl',
            marker: {
                size: state.markerSize,
                color: expr,
                colorscale: EXPRESSION_COLORSCALE,
                opacity: state.opacity,
                colorbar: {
                    title: state.selectedGene,
                    thickness: 15,
                    len: 0.5,
                }
            },
            hovertemplate: `x: %{x}<br>y: %{y}<br>${state.selectedGene}: %{marker.color:.2f}<extra></extra>`
        }];

        updateLegendForExpression();

    } else {
        // Categorical mode
        const categories = [...new Set(colorByData)].sort();
        const colorMap = {};
        categories.forEach((cat, i) => {
            colorMap[cat] = CATEGORICAL_COLORS[i % CATEGORICAL_COLORS.length];
        });

        const colors = colorByData.map(cat => colorMap[cat]);

        // Apply highlighting filter
        let opacity = colorByData.map(() => state.opacity);
        if (state.highlightedCategories.size > 0) {
            opacity = colorByData.map(cat =>
                state.highlightedCategories.has(cat) ? state.opacity : 0.05
            );
        }

        data = [{
            x: coordinates.x,
            y: coordinates.y,
            mode: 'markers',
            type: 'scattergl',
            marker: {
                size: state.markerSize,
                color: colors,
                opacity: opacity,
            },
            text: colorByData,
            hovertemplate: `x: %{x}<br>y: %{y}<br>${state.colorBy}: %{text}<extra></extra>`
        }];

        updateLegendForCategories(categories, colorMap, colorByData);
    }

    layout = {
        xaxis: {
            title: 'X',
            scaleanchor: 'y',
            scaleratio: 1,
        },
        yaxis: {
            title: 'Y',
        },
        margin: { t: 20, b: 40, l: 50, r: 20 },
        hovermode: 'closest',
        dragmode: 'pan',
    };

    const config = {
        responsive: true,
        displayModeBar: true,
        modeBarButtonsToRemove: ['select2d', 'lasso2d'],
        displaylogo: false,
    };

    Plotly.react(elements.plot, data, layout, config);
}

function updateLegendForCategories(categories, colorMap, colorByData) {
    // Count cells per category
    const counts = {};
    colorByData.forEach(cat => {
        counts[cat] = (counts[cat] || 0) + 1;
    });

    elements.legendTitle.textContent = formatColumnName(state.colorBy);

    elements.legendList.innerHTML = categories.map(cat => {
        const isDimmed = state.highlightedCategories.size > 0 &&
            !state.highlightedCategories.has(cat);
        return `
            <div class="legend-item ${isDimmed ? 'dimmed' : ''}" data-category="${cat}">
                <div class="legend-color" style="background: ${colorMap[cat]}"></div>
                <span class="legend-label" title="${cat}">${cat}</span>
                <span class="legend-count">${counts[cat]?.toLocaleString() || 0}</span>
            </div>
        `;
    }).join('');

    // Add click handlers for legend items
    elements.legendList.querySelectorAll('.legend-item').forEach(item => {
        item.addEventListener('click', () => {
            const cat = item.dataset.category;
            toggleCategoryHighlight(cat);
        });
    });
}

function updateLegendForExpression() {
    elements.legendTitle.textContent = `Gene: ${state.selectedGene}`;

    elements.legendList.innerHTML = `
        <div style="padding: 20px; text-align: center; color: #666;">
            <div style="margin-bottom: 10px;">Expression Level</div>
            <div style="height: 20px; border-radius: 4px; background: linear-gradient(to right, #e0e0e0, #ffeda0, #feb24c, #f03b20, #bd0026);"></div>
            <div style="display: flex; justify-content: space-between; font-size: 0.8rem; margin-top: 5px;">
                <span>Low</span>
                <span>High</span>
            </div>
        </div>
    `;
}

function toggleCategoryHighlight(category) {
    if (state.highlightedCategories.has(category)) {
        state.highlightedCategories.delete(category);
    } else {
        state.highlightedCategories.add(category);
    }

    // If all are deselected, show all
    const allCategories = [...new Set(state.metadata.metadata[state.colorBy])];
    if (state.highlightedCategories.size === allCategories.length) {
        state.highlightedCategories.clear();
    }

    renderPlot();
}

function setupEventListeners() {
    // Dataset change
    elements.datasetSelect.addEventListener('change', (e) => {
        loadDataset(e.target.value);
    });

    // Color by change
    elements.colorBySelect.addEventListener('change', (e) => {
        state.colorBy = e.target.value;
        state.highlightedCategories.clear();
        renderPlot();
    });

    // Marker size
    elements.markerSize.addEventListener('input', (e) => {
        state.markerSize = parseInt(e.target.value);
        elements.markerSizeValue.textContent = state.markerSize;
        renderPlot();
    });

    // Opacity
    elements.opacity.addEventListener('input', (e) => {
        state.opacity = parseInt(e.target.value) / 100;
        elements.opacityValue.textContent = `${e.target.value}%`;
        renderPlot();
    });

    // Gene search
    elements.geneSearch.addEventListener('input', debounce((e) => {
        const query = e.target.value.trim().toUpperCase();
        if (query.length < 2) {
            hideAutocomplete();
            return;
        }

        const genes = state.metadata.top_genes || [];
        const matches = genes.filter(g => g.toUpperCase().includes(query)).slice(0, 20);
        showAutocomplete(matches, query);
    }, 200));

    elements.geneSearch.addEventListener('blur', () => {
        setTimeout(hideAutocomplete, 200);
    });

    elements.geneSearch.addEventListener('keydown', (e) => {
        if (e.key === 'Escape') {
            hideAutocomplete();
            if (!state.selectedGene) {
                e.target.value = '';
            }
        } else if (e.key === 'Enter') {
            const first = elements.geneAutocomplete.querySelector('.autocomplete-item');
            if (first) {
                selectGene(first.dataset.gene);
            }
        }
    });

    // Reset zoom
    elements.btnResetZoom.addEventListener('click', () => {
        Plotly.relayout(elements.plot, {
            'xaxis.autorange': true,
            'yaxis.autorange': true,
        });
    });

    // Download
    elements.btnDownload.addEventListener('click', () => {
        Plotly.downloadImage(elements.plot, {
            format: 'png',
            width: 1600,
            height: 1200,
            filename: `${state.currentDataset}_${state.colorBy}`
        });
    });
}

function showAutocomplete(matches, query) {
    if (matches.length === 0) {
        hideAutocomplete();
        return;
    }

    elements.geneAutocomplete.innerHTML = matches.map(gene => {
        const highlighted = gene.replace(
            new RegExp(`(${query})`, 'gi'),
            '<mark>$1</mark>'
        );
        return `<div class="autocomplete-item" data-gene="${gene}">${highlighted}</div>`;
    }).join('');

    elements.geneAutocomplete.classList.add('show');

    // Add click handlers
    elements.geneAutocomplete.querySelectorAll('.autocomplete-item').forEach(item => {
        item.addEventListener('click', () => {
            selectGene(item.dataset.gene);
        });
    });
}

function hideAutocomplete() {
    elements.geneAutocomplete.classList.remove('show');
}

function selectGene(gene) {
    state.selectedGene = gene;
    elements.geneSearch.value = gene;
    hideAutocomplete();
    state.highlightedCategories.clear();
    renderPlot();
}

function showLoading(show, text = 'Loading data...') {
    elements.loading.style.display = show ? 'flex' : 'none';
    elements.loading.querySelector('.loading-text').textContent = text;
}

function showError(message) {
    elements.legendList.innerHTML = `
        <div style="padding: 20px; color: #e74c3c; text-align: center;">
            ⚠️ ${message}
        </div>
    `;
}

function debounce(fn, delay) {
    let timeout;
    return (...args) => {
        clearTimeout(timeout);
        timeout = setTimeout(() => fn(...args), delay);
    };
}

// Start
document.addEventListener('DOMContentLoaded', init);
