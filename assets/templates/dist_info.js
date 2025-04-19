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
    const currentTotal = shadow + sun;

    // Calcular la distancia máxima entre todas las rutas
    const maxDistance = Math.max(...routeDistances.map(r => r.distance_shadow + r.distance_sun));
    
    // El porcentaje total ocupado por la ruta actual respecto a la ruta más larga
    const filledPercentage = maxDistance ? (currentTotal / maxDistance) * 100 : 0;
    
    // Proporción interna para sombra y sol
    const shadowPct = currentTotal ? (shadow / currentTotal) * filledPercentage : 0;
    const sunPct = currentTotal ? (sun / currentTotal) * filledPercentage : 0;
    
    // Porcentaje restante en gris
    const remainingPct = 100 - filledPercentage;

    const shadowSegment = document.getElementById('shadowSegment');
    const sunSegment = document.getElementById('sunSegment');
    const remainingSegment = document.getElementById('remainingSegment');

    shadowSegment.style.width = shadowPct + '%';
    shadowSegment.style.background = '#1795d4';
    shadowSegment.textContent = shadow;

    sunSegment.style.width = sunPct + '%';
    sunSegment.style.background = '#ff271c';
    sunSegment.textContent = sun;
    
    remainingSegment.style.width = remainingPct + '%';
    // Mantener fondo gris para el segmento restante
    remainingSegment.textContent = ''; // O puedes mostrar el valor de diferencia si lo prefieres

    infoDiv.style.display = 'block';
}

function hideMeterInfo() {
    const infoDiv = document.getElementById('meterInfo');
    if (infoDiv) {
        infoDiv.style.display = 'none';
    }
}