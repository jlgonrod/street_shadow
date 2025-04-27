document.addEventListener('DOMContentLoaded', function() {
    const slider = document.getElementById('routeSlider');
    // If there is only one route (max ≤ 1), replace the slider and the "Show All Routes" checkbox
    if (slider && parseInt(slider.getAttribute('max')) <= 1) {
        const sliderContainer = document.getElementById('sliderContainer');
        sliderContainer.innerHTML = `
            <div style="text-align: center;">Only one route is available,<br>selection is not possible.</div>
            <br>
            <input type="checkbox" id="toggleShadow" checked onchange="updateRoutes()">
            <label for="toggleShadow">Show Shadows</label>
        `;
        // Update routes, forcing showAllRoutes to false
        updateRoutes();
    }
});

function updateRoutes() {
    const slider = document.getElementById('routeSlider');
    // If the "showAllRoutes" checkbox does not exist, simulate it as unchecked (false)
    let showAll = document.getElementById('showAllRoutes');
    if (!showAll) {
        showAll = { checked: false };
    }
    const toggleShadow = document.getElementById('toggleShadow');
    // If the slider does not exist (single route), use 0 by default
    const selectedRoute = slider ? slider.value : 0;
    // Get all path elements from the map
    const paths = document.querySelectorAll('path');
    let routeIndex = 0;
    paths.forEach(path => {
        const stroke = path.getAttribute('stroke');
        // If the path is a route (stroke is not black)
        if (stroke && stroke.toLowerCase() !== '#000000' && stroke.toLowerCase() !== 'black') {
            if (showAll.checked) {
                path.style.display = '';
            } else {
                path.style.display = (routeIndex == selectedRoute) ? '' : 'none';
            }
            routeIndex++;
        } else {
            // For shadows (black stroke), check the state of the shadow checkbox
            path.style.display = toggleShadow.checked ? '' : 'none';
        }
    });
    if (!showAll.checked) {
        // Force the info container to be displayed so its offsetWidth is calculated correctly
        const infoDiv = document.getElementById('timeInfo');
        if (infoDiv) {
            infoDiv.style.display = 'block';
        }
        // Update the information in the next frame to ensure the layout is updated
        window.requestAnimationFrame(() => {
            updateTimeInfo(selectedRoute, showAll.checked);
        });
    } else {
        hideTimeInfo();
    }
}

function sliderChanged() {
    // When the slider changes, deactivate "Show All Routes" if it exists
    const showAllEl = document.getElementById('showAllRoutes');
    if (showAllEl) {
       showAllEl.checked = false;
    }
    updateRoutes();
}