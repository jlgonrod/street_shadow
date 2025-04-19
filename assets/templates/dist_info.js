function updateMeterInfo(selectedIndex) {
    if (typeof routeDistances === 'undefined') return;
    const infoDiv = document.getElementById('meterInfo');
    if (!infoDiv) return;

    const routeData = routeDistances[selectedIndex];
    if (!routeData) {
        infoDiv.style.display = 'none';
        return;
    }
    const shadow = routeData.distance_shadow;
    const sun = routeData.distance_sun;
    const total = shadow + sun;
    // Evitar división por cero
    const shadowPerc = total ? (shadow / total) * 100 : 0;
    const sunPerc = total ? (sun / total) * 100 : 0;

    const shadowSegment = document.getElementById('shadowSegment');
    const sunSegment = document.getElementById('sunSegment');
    shadowSegment.style.width = shadowPerc + '%';
    sunSegment.style.width = sunPerc + '%';

    // Mostrar los valores numéricos dentro de cada segmento
    shadowSegment.textContent = shadow;
    sunSegment.textContent = sun;

    infoDiv.style.display = 'block';
}

function hideMeterInfo() {
    const infoDiv = document.getElementById('meterInfo');
    if (infoDiv) {
        infoDiv.style.display = 'none';
    }
}