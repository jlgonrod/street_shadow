document.addEventListener('DOMContentLoaded', function() {
    const slider = document.getElementById('routeSlider');
    // Si solo hay una ruta (max ≤ 1) se reemplaza el slider y el checkbox "Show All Routes"
    if (slider && parseInt(slider.getAttribute('max')) <= 1) {
        const sliderContainer = document.getElementById('sliderContainer');
        sliderContainer.innerHTML = `
            <div style="text-align: center;">Only one route is available,<br>selection is not possible.</div>
            <br>
            <input type="checkbox" id="toggleShadow" checked onchange="updateRoutes()">
            <label for="toggleShadow">Show Shadows</label>
        `;
        // Actualizamos las rutas, forzando showAllRoutes como false
        updateRoutes();
    }
});

function updateRoutes() {
    const slider = document.getElementById('routeSlider');
    // Si no existe el checkbox "showAllRoutes", se simula que está desactivado (false)
    let showAll = document.getElementById('showAllRoutes');
    if (!showAll) {
        showAll = { checked: false };
    }
    const toggleShadow = document.getElementById('toggleShadow');
    // Si el slider no existe (ruta única), se usa 0 por defecto
    const selectedRoute = slider ? slider.value : 0;
    // Obtiene todos los elementos path del mapa
    const paths = document.querySelectorAll('path');
    let routeIndex = 0;
    paths.forEach(path => {
        const stroke = path.getAttribute('stroke');
        // Si el path es una ruta (stroke distinto de negro)
        if (stroke && stroke.toLowerCase() !== '#000000' && stroke.toLowerCase() !== 'black') {
            if (showAll.checked) {
                path.style.display = '';
            } else {
                path.style.display = (routeIndex == selectedRoute) ? '' : 'none';
            }
            routeIndex++;
        } else {
            // Para sombras (stroke negro), se verifica el estado del checkbox de sombras
            path.style.display = toggleShadow.checked ? '' : 'none';
        }
    });
    if (!showAll.checked) {
        // Forzamos que se muestre el contenedor de información para que se calcule el offsetWidth correctamente
        const infoDiv = document.getElementById('timeInfo');
        if (infoDiv) {
            infoDiv.style.display = 'block';
        }
        // Actualiza la información en el siguiente frame para garantizar que el layout esté actualizado
        window.requestAnimationFrame(() => {
            updateTimeInfo(selectedRoute, showAll.checked);
        });
    } else {
        hideTimeInfo();
    }
}

function sliderChanged() {
    // Al cambiar el slider se desactiva "Show All Routes", si este existe
    const showAllEl = document.getElementById('showAllRoutes');
    if (showAllEl) {
       showAllEl.checked = false;
    }
    updateRoutes();
}