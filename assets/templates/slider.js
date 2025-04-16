function updateRoutes() {
    const slider = document.getElementById('routeSlider');
    const showAll = document.getElementById('showAllRoutes');
    const toggleShadow = document.getElementById('toggleShadow');
    const selectedRoute = slider.value;
    // Get all SVG path elements on the map
    const paths = document.querySelectorAll('path');
    let routeIndex = 0;
    paths.forEach(path => {
        const stroke = path.getAttribute('stroke');
        // if the path is a route (non-black stroke)
        if (stroke && stroke.toLowerCase() !== '#000000' && stroke.toLowerCase() !== 'black') {
            if (showAll.checked) {
                path.style.display = '';
            } else {
                path.style.display = (routeIndex == selectedRoute) ? '' : 'none';
            }
            routeIndex++;
        } else {
            // for shadows (drawn in black), check the toggleShadow status
            if (toggleShadow.checked) {
                path.style.display = '';
            } else {
                path.style.display = 'none';
            }
        }
    });
}

function sliderChanged() {
    // Automatically uncheck "Show All Routes" when slider changes
    document.getElementById('showAllRoutes').checked = false;
    updateRoutes();
}