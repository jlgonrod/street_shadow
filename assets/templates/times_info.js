function updateTimeInfo(selectedIndex, showAllRoutes) {
    if (typeof routeTimes === 'undefined') return;
    const infoDiv = document.getElementById('timeInfo');
    if (!infoDiv) return;

    const routeData = routeTimes[selectedIndex];
    if (!routeData) {
        infoDiv.style.display = 'none';
        return;
    }

    // Valores de la ruta seleccionada
    const shadow = routeData.time_shadow;
    const sun = routeData.time_sun;
    const total = shadow + sun;

    let maxShadow, maxSun;
    if (showAllRoutes) {
        // Usamos máximos globales
        maxShadow = Math.max(...routeTimes.map(r => r.time_shadow));
        maxSun = Math.max(...routeTimes.map(r => r.time_sun));
    } else {
        // Usamos los valores de la ruta seleccionada
        maxShadow = shadow;
        maxSun = sun;
    }

    // Usar el ancho total del contenedor (la barra gris) como referencia
    const timeBar = document.getElementById('timeBar');
    const totalWidth = timeBar.offsetWidth;

    // Calcular el factor de escala basado en los máximos globales o seleccionados
    const scale = totalWidth / (maxShadow + maxSun);

    // El origen se ubica en el centro de la barra
    const origin = (maxShadow / (maxShadow + maxSun)) * totalWidth;

    const shadowSegment = document.getElementById('timeShadowSegment');
    const sunSegment = document.getElementById('timeSunSegment');
    const divider = document.getElementById('timeDivider');

    // Configurar el segmento de Sombra
    const shadowWidth = shadow * scale;
    shadowSegment.style.width = shadowWidth + 'px';
    // Ubicamos su extremo derecho en el origen
    shadowSegment.style.left = (origin - shadowWidth) + 'px';
    shadowSegment.textContent = shadow + ' min';

    // Configurar el segmento de Sol
    const sunWidth = sun * scale;
    sunSegment.style.width = sunWidth + 'px';
    sunSegment.style.left = origin + 'px';
    sunSegment.textContent = sun + ' min';

    // Posicionar el divisor en el origen
    divider.style.left = origin + 'px';
    divider.style.transform = 'none';

    // Mostrar el total en la parte inferior
    const totalTimeDiv = document.getElementById('totalTime');
    totalTimeDiv.textContent = 'Total min: ' + total;

    infoDiv.style.display = 'block';
}

function hideTimeInfo() {
    const infoDiv = document.getElementById('timeInfo');
    if (infoDiv) {
        infoDiv.style.display = 'none';
    }
}