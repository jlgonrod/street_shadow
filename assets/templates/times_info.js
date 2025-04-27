function updateTimeInfo(selectedIndex, showAllRoutes) {
    if (typeof routeTimes === 'undefined') return;
    const infoDiv = document.getElementById('timeInfo');
    if (!infoDiv) return;

    const routeData = routeTimes[selectedIndex];
    if (!routeData) {
        infoDiv.style.display = 'none';
        return;
    }

    // Values of the selected route
    const shadow = routeData.time_shadow;
    const sun = routeData.time_sun;
    const total = shadow + sun;

    // Always use global maximums to keep bars the same size
    const maxShadow = Math.max(...routeTimes.map(r => r.time_shadow));
    const maxSun = Math.max(...routeTimes.map(r => r.time_sun));

    // Use the total width of the container (the gray bar) as a reference
    const timeBar = document.getElementById('timeBar');
    const totalWidth = timeBar.offsetWidth;

    // Calculate the scaling factor using the global maximums
    const scale = totalWidth / (maxShadow + maxSun);

    // The origin is located at the center of the bar based on the maximums
    const origin = (maxShadow / (maxShadow + maxSun)) * totalWidth;

    const shadowSegment = document.getElementById('timeShadowSegment');
    const sunSegment = document.getElementById('timeSunSegment');
    const divider = document.getElementById('timeDivider');

    // Configure the Shadow segment
    const shadowWidth = shadow * scale;
    shadowSegment.style.width = shadowWidth + 'px';
    // Position its right end at the origin
    shadowSegment.style.left = (origin - shadowWidth) + 'px';
    shadowSegment.textContent = shadow ? shadow + ' min' : '';

    // Configure the Sun segment
    const sunWidth = sun * scale;
    sunSegment.style.width = sunWidth + 'px';
    sunSegment.style.left = origin + 'px';
    sunSegment.textContent = sun ? sun + ' min' : '';

    // Position the divider at the origin
    divider.style.left = origin + 'px';
    divider.style.transform = 'none';

    // Display the total at the bottom
    const totalTimeDiv = document.getElementById('totalTime');
    totalTimeDiv.textContent = total ? 'Total min: ' + total : '';

    infoDiv.style.display = 'block';
}

function hideTimeInfo() {
    const infoDiv = document.getElementById('timeInfo');
    if (infoDiv) {
        infoDiv.style.display = 'none';
    }
}