// ------------------------ ------------------------
// BASIC SETUP
// ------------------------------------------------
// Create an empty scene
var scene = new THREE.Scene();

// Create a basic perspective camera
var camera = new THREE.PerspectiveCamera( 75, window.innerWidth/window.innerHeight, 0.1, 10000000 );
camera.position.z = 500;

// Create a renderer with antialiasing
var renderer = new THREE.WebGLRenderer({antialias:true});

// White background
renderer.setClearColor("#ffffff");

// render to entire canvas
renderer.setSize( window.innerWidth, window.innerHeight );

// enable keyboard event handler:
window.addEventListener( "keydown", keyboardHandler )

// populate the info overlay once the page loads:
window.addEventListener( "DOMContentLoaded", update_overlay )

// Append Renderer to DOM
document.body.appendChild( renderer.domElement );

// timekeeping for animation
var clock = new THREE.Clock();
var t_elapsed = 0


/////////////////
// Scene Setup //
/////////////////

// "constants"
var SCALE_FUDGE = 0.01 // scale the scene by this much
var EIG_T_FUDGE = 100 // scale the translation perturbations by this much
var EIG_R_FUDGE = 100 // scale the rotation perturbations by this much
var CURRENT_EIG = 0 // the curent eigenvector
var CAM_PT_SCALE = 1.0
var DEF_SCALE_PCTILE = 0.5

var n_cameras = sample_scene["num_cameras"]
var n_points = sample_scene["num_points"]
var n_eigs = sample_scene["eigenvectors"].length

var positions = []; // position parameters of all cameras
var rotations = []; // rotation parameters of all cameras

// array to store the perturbations of each camera parameter
// from the current eigenvector:
var wiggles = Array(n_cameras*6).fill(0);

// set default scale to zero (will be auto-calculated) for any cameras
// it's not set for in the html file
while (DEFAULT_SCALES.length < n_cameras) {
    DEFAULT_SCALES.push(0)
}


// set up initial camera positions and rotations
for ( var i = 0; i < n_cameras; i ++ ) {
    cp = sample_scene["camera_parameters"];
    cc = get_pos( cp[9*i], cp[9*i+1], cp[9*i+2], cp[9*i+3], cp[9*i+4], cp[9*i+5] );

    positions.push( cc.x * SCALE_FUDGE );
    positions.push( cc.y * SCALE_FUDGE );
    positions.push( cc.z * SCALE_FUDGE );

    rotations.push( cp[9*i+3] );
    rotations.push( cp[9*i+4] );
    rotations.push( cp[9*i+5] );
}

//var axesHelper = new THREE.AxesHelper( 5 );
//scene.add( axesHelper );

/* Create a frustum for each camera */
frusta = []

for (var i = 0; i < n_cameras; i++ ) {
    var wire_material = new THREE.LineBasicMaterial( { color: 0x000000, linewidth: 1 } );
    var cone_geom = new THREE.EdgesGeometry( new THREE.ConeGeometry(1, 2, 4, 1, false ) )
    var wireframe = new THREE.LineSegments(cone_geom, wire_material )
    // set frustum to the canonical orientation looking down the -z axis:
    wireframe.rotateX(Math.PI/2);
    wireframe.rotateY(Math.PI/4);

    // with its pointy end at the origin:
    wireframe.translateY(-1);

    // wrap the above in a wrapper object that can be translated and rotated to
    // place it in the scene:
    var frustum = new THREE.Object3D();
    frustum.add(wireframe);
    frusta.push(frustum);
    scene.add(frustum)

    // position and orient the camera
    frustum.position.set(positions[3*i], positions[3*i+1], positions[3*i+2])

    var r1 = rotations[3*i]
    var r2 = rotations[3*i+1]
    var r3 = rotations[3*i+2]

    axis = get_axis(r1, r2, r3);
    angle = get_angle(r1, r2, r3);
    frustum.setRotationFromAxisAngle(axis, angle);
}


// set up scene point positions
var pt_positions = []
for ( var i = 0; i < n_points; i ++ ) {
    pt_positions.push(SCALE_FUDGE * sample_scene["point_parameters"][3*i]);
    pt_positions.push(SCALE_FUDGE * sample_scene["point_parameters"][3*i+1]);
    pt_positions.push(SCALE_FUDGE * sample_scene["point_parameters"][3*i+2]);
}

set_eigenvector( CURRENT_EIG );

pts_geom = new THREE.BufferGeometry();
pts_geom.addAttribute( 'position', new THREE.Float32BufferAttribute( pt_positions, 3 ) );

var scene_points = new THREE.Points( pts_geom, new THREE.PointsMaterial( { size: 0.5, color: 0x226688 } ) );
scene.add( scene_points );

var controls = new THREE.OrbitControls( camera );
controls.screenSpacePanning = true;



/** Set the eigenvector to eig_index. Sets up the wiggles array
 * to contain the component of the eigenvector for each parameter,
 * which is later used in wiggle_cameras to vary the camera locations
 * and orientations. */
function set_eigenvector( eig_index ) {
    CURRENT_EV = sample_scene["eigenvalues"][eig_index]
    var t_mags = []

    for ( var i = 0; i < n_cameras; i++ ) {
        er1 = sample_scene["eigenvectors"][eig_index][6*i+0]
        er2 = sample_scene["eigenvectors"][eig_index][6*i+1]
        er3 = sample_scene["eigenvectors"][eig_index][6*i+2]
        et1 = sample_scene["eigenvectors"][eig_index][6*i+3]
        et2 = sample_scene["eigenvectors"][eig_index][6*i+4]
        et3 = sample_scene["eigenvectors"][eig_index][6*i+5]

        wiggles[6*i+0] = er1;
        wiggles[6*i+1] = er2;
        wiggles[6*i+2] = er3;

        wiggles[6*i+3] = et1;
        wiggles[6*i+4] = et2;
        wiggles[6*i+5] = et3;

        t_mags.push(Math.sqrt(et1*et1 + et2*et2 + et3*et3));
    }

    if ( DEFAULT_SCALES[eig_index] == 0 ) {
        t_mags.sort()

        pctile = t_mags[Math.floor(DEF_SCALE_PCTILE * n_cameras)]
        EIG_T_FUDGE = pctile * 1e8
        EIG_R_FUDGE = pctile * 1e8
        DEFAULT_SCALES[eig_index] = pctile * 1e8 // memorize for later
    } else {
        EIG_T_FUDGE = DEFAULT_SCALES[eig_index]
        EIG_R_FUDGE = DEFAULT_SCALES[eig_index]
    }


}


/** Get the unit-vector axis from an axis-angle vector */
function get_axis(r1, r2, r3) {
    rnorm = Math.sqrt(r1 * r1 + r2 * r2 + r3 * r3)
    k1 = r1 / rnorm
    k2 = r2 / rnorm
    k3 = r3 / rnorm
    return new THREE.Vector3(k1, k2, k3)
}

/** Get the angle from an axis-angle vector */
function get_angle(r1, r2, r3) {
    return Math.sqrt(r1 * r1 + r2 * r2 + r3 * r3)
}

/* Get position of a camera center, which is -R'*t
 * Takes the axis-angle vector [r1 r2 r3]' and
 * translation vector [t1 t2 t3]' */
function get_pos( r1, r2, r3, t1, t2, t3 ) {
    rnorm = Math.sqrt(r1 * r1 + r2 * r2 + r3 * r3)
    k1 = r1 / rnorm
    k2 = r2 / rnorm
    k3 = r3 / rnorm

    k = new THREE.Vector3( k1, k2, k3 )
    t = new THREE.Vector3( t1, t2, t3 )
    t.applyAxisAngle( k, -rnorm )
    t.multiplyScalar(-1)

    return t
}

/** Implement keyboard interaction */
function keyboardHandler( e ) {
    if ( e.key == "h" ) {
        EIG_T_FUDGE *= 2;
    } else if ( e.key == "l" ) {
        EIG_T_FUDGE /= 2;
    } else if ( e.key == "g" ) {
        EIG_R_FUDGE *= 2;
    } else if ( e.key == ";" ) {
        EIG_R_FUDGE /= 2;
    } else if ( e.key == "j") {
        CURRENT_EIG = (CURRENT_EIG + 1) % n_eigs;
        set_eigenvector(CURRENT_EIG);
    } else if ( e.key == "k") {
        CURRENT_EIG -= 1;
        if ( CURRENT_EIG < 0 ) {
            CURRENT_EIG = n_eigs - 1;
        }
        set_eigenvector(CURRENT_EIG);
    } else if ( e.key == 'p' ) {
        scene_points.visible = !scene_points.visible; // toggle scene point visibility
    } else if ( e.key == '?' ) {
        var overlay = document.getElementById("kb_overlay");
        if (overlay.style.visibility == "hidden") {
            overlay.style.visibility = "visible";
        } else {
            overlay.style.visibility = "hidden";
        }
    }
    update_overlay();
}

/** Update the text overlay showing current display setting */
function update_overlay() {
    var overlay = document.getElementById("overlay")
    display_html = "<table>";
    display_html += table_row("Eigenvector:", CURRENT_EIG);
    display_html += table_row("Eigenvalue:", CURRENT_EV);
    display_html += table_row("Eig T scale:", EIG_T_FUDGE);
    display_html += table_row("Eig R scale:", EIG_R_FUDGE);
    display_html += "</table>";
    overlay.innerHTML = display_html;
}

function table_row(label, value) {
    return "<tr><td>" + label + "</td><td>" + value + "</td></tr>";
}

/** Update function to wiggle points by their wiggle vectors * sin(t) */
function wiggle_cameras(t) {

    sint = Math.sin(t); // current wiggle amount

    for ( var i = 0; i < n_cameras; i++ ) {

        var cam = frusta[i];

        var x = positions[3*i] + EIG_T_FUDGE * wiggles[6*i+3] * sint;
        var y = positions[3*i+1] + EIG_T_FUDGE * wiggles[6*i+4] * sint;
        var z = positions[3*i+2] + EIG_T_FUDGE * wiggles[6*i+5] * sint;

        var r1 = rotations[3*i] + EIG_R_FUDGE * wiggles[6*i] * sint;
        var r2 = rotations[3*i+1] + EIG_R_FUDGE * wiggles[6*i+1] * sint;
        var r3 = rotations[3*i+2] + EIG_R_FUDGE * wiggles[6*i+2] * sint;

        axis = get_axis(r1, r2, r3);
        angle = get_angle(r1, r2, r3);

        cam.position.set(x, y, z);
        cam.setRotationFromAxisAngle(axis, angle)
    }
}

/////////////////
// Render Loop //
/////////////////
var render = function () {
    requestAnimationFrame( render );

    t_delta = clock.getDelta()
    t_elapsed = clock.getElapsedTime()

    controls.update( t_delta )

    // Render the scene
    renderer.render( scene, camera );

    wiggle_cameras( t_elapsed * 4 );

};

render();
