function updateTimeInfo(selectedIndex) {
    if (typeof routeTimes === 'undefined') return;
    const infoDiv = document.getElementById('timeInfo');
    if (!infoDiv) return;

    const routeData = routeTimes[selectedIndex];
    if (!routeData) {
        infoDiv.style.display = 'none';
        return; // <---- Se agrega el return para detener la ejecución
    }
    
    const shadow = routeData.time_shadow;
    const sun = routeData.time_sun;
    const currentTotal = shadow + sun;

    // Calcular el tiempo máximo entre todas las rutas
    const maxTime = Math.max(...routeTimes.map(r => r.time_shadow + r.time_sun));
    
    // Porcentaje que ocupa la ruta actual respecto a la ruta más larga
    const filledPercentage = maxTime ? (currentTotal / maxTime) * 100 : 0;
    
    // Proporciones para los segmentos de sombra y sol
    const shadowPct = currentTotal ? (shadow / currentTotal) * filledPercentage : 0;
    const sunPct = currentTotal ? (sun / currentTotal) * filledPercentage : 0;
    
    // Porcentaje restante para completar 100%
    const remainingPct = 100 - filledPercentage;

    const shadowSegment = document.getElementById('timeShadowSegment');
    const sunSegment = document.getElementById('timeSunSegment');
    const remainingSegment = document.getElementById('timeRemainingSegment');

    shadowSegment.style.width = shadowPct + '%';
    shadowSegment.style.background = '#1795d4';
    shadowSegment.textContent = shadow;

    sunSegment.style.width = sunPct + '%';
    sunSegment.style.background = '#ff271c';
    sunSegment.textContent = sun;
    
    remainingSegment.style.width = remainingPct + '%';
    remainingSegment.style.background = '#cccccc';
    remainingSegment.textContent = '';
    
    infoDiv.style.display = 'block';
}

function hideTimeInfo() {
    const infoDiv = document.getElementById('timeInfo');
    if (infoDiv) {
        infoDiv.style.display = 'none';
    }
}