// ------------------------------------------------
// BASIC SETUP
// ------------------------------------------------
//var THREE = require('three')
//import { FirstPersonControls } from 'three/examples/jsm/controls/FirstPersonControls.js';
//import * as THREE from 'three.module.js';
//
// Create an empty scene
var scene = new THREE.Scene();

// Create a basic perspective camera
var camera = new THREE.PerspectiveCamera( 75, window.innerWidth/window.innerHeight, 0.1, 10000000 );
camera.position.z = 500;

// Create a renderer with Antialiasing
var renderer = new THREE.WebGLRenderer({antialias:true});

// Configure renderer clear color
renderer.setClearColor("#ffffff");

// Configure renderer size
renderer.setSize( window.innerWidth, window.innerHeight );

window.addEventListener( "keydown", keyboardHandler )
window.addEventListener( "DOMContentLoaded", update_overlay )


// Append Renderer to DOM
document.body.appendChild( renderer.domElement );


var clock = new THREE.Clock();
var t_elapsed = 0


var capturer = new CCapture( {
    format: 'png',
    framerate: 120,
    verbose: false,
    timeLimit: 1,
    autoSaveTime: 10,
    name: get_timestamp()
} );

/////////////////
// Scene Setup //
/////////////////

// "constants"
var SCALE_FUDGE = 0.01
var EIG_T_FUDGE = 100
var EIG_R_FUDGE = 100
var CURRENT_EIG = 0
var STATIC = true
var CAM_PT_SCALE = 1.0

var n_cameras = sample_scene["num_cameras"]
var n_points= sample_scene["num_points"]
var n_eigs = sample_scene["eigenvectors"].length

console.log("Cameras: " + n_cameras)
console.log("Points: " + n_points)
console.log("Eigenvectors: " + n_eigs)
//var cameras = 10
var positions = [];
var rotations = []
// var colors = [];
// var sizes = [];
var wiggles = Array(n_cameras*6).fill(0);
//var color = new THREE.Color();



// set up initial camera positions
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

frusta = []

for (var i = 0; i < n_cameras; i++ ) {
    var wire_material = new THREE.LineBasicMaterial( { color: 0x000000, linewidth: 1 } );
    var cone_geom = new THREE.EdgesGeometry( new THREE.ConeGeometry(1, 2, 4, 1, false ) )
    var wireframe = new THREE.LineSegments(cone_geom, wire_material )
    wireframe.rotateX(Math.PI/2);
    wireframe.rotateY(Math.PI/4);
    wireframe.translateY(-1);
    var frustum = new THREE.Object3D();
    frustum.add(wireframe)
    //var pos = new THREE.Vector3();
    //pos.fromArray(positions, 3*i)
    frusta.push(frustum);
    scene.add(frustum)

    frustum.position.set(positions[3*i], positions[3*i+1], positions[3*i+2])

    var r1 = rotations[3*i]
    var r2 = rotations[3*i+1]
    var r3 = rotations[3*i+2]

    axis = get_axis(r1, r2, r3);
    angle = get_angle(r1, r2, r3);
    frustum.setRotationFromAxisAngle(axis, angle);
    //var R = new THREE.Matrix4
    //R.makeRotationAxis(axis, angle)
    //wireframe.matrix.copy(R)
    //console.log(wireframe.matrix)
    //console.log(wireframe.position);
    //console.log(wireframe.rotation);
}


// set up scene point positions
var pt_positions = []
for ( var i = 0; i < n_points; i ++ ) {
    
    pt_positions.push(SCALE_FUDGE * sample_scene["point_parameters"][3*i])
    pt_positions.push(SCALE_FUDGE * sample_scene["point_parameters"][3*i+1])
    pt_positions.push(SCALE_FUDGE * sample_scene["point_parameters"][3*i+2])
}

set_eigenvector( CURRENT_EIG );


// set up geometry and create the Points object
// color and size seem to be ignored, not sure what's going on there
geometry = new THREE.BufferGeometry();
geometry.addAttribute( 'position', new THREE.Float32BufferAttribute( positions, 3 ) );
//geometry.addAttribute( 'color', new THREE.Float32BufferAttribute( colors, 3 ) );
//geometry.addAttribute( 'size', new THREE.Float32BufferAttribute( sizes, 1 ) );
geometry.attributes.position.needsUpdate = true;

pts_geom = new THREE.BufferGeometry();
pts_geom.addAttribute( 'position', new THREE.Float32BufferAttribute( pt_positions, 3 ) );


var camera_points = new THREE.Points( geometry, new THREE.PointsMaterial( { size: 0.5, color: 0x000000} ) );
var scene_points = new THREE.Points( pts_geom, new THREE.PointsMaterial( { size: 0.5, color: 0x226688 } ) );
scene.add( camera_points );
scene.add( scene_points );


var controls = new THREE.OrbitControls( camera );
controls.screenSpacePanning = true;

function set_eigenvector( eig_index ) {
    //max_t = 0
    //max_r = 0
    CURRENT_EV = sample_scene["eigenvalues"][eig_index]
    for ( var i = 0; i < n_cameras; i++ ) {
        er1 = sample_scene["eigenvectors"][eig_index][6*i+0]
        er2 = sample_scene["eigenvectors"][eig_index][6*i+1]
        er3 = sample_scene["eigenvectors"][eig_index][6*i+2]
        et1 = sample_scene["eigenvectors"][eig_index][6*i+3]
        et2 = sample_scene["eigenvectors"][eig_index][6*i+4]
        et3 = sample_scene["eigenvectors"][eig_index][6*i+5]

        //max_t = Math.max(max_t, et1, et2, et3)
        //max_r = Math.max(max_r, er1, er2, er3)

        //if (max_t == et1 || max_t == et2 || max_t == et3) {
            //max_ti = i
        //}
        //if (max_r == er1 || max_r == er2 || max_r == er3) {
            //max_ri = i
        //}

        wiggles[6*i+0] = er1;
        wiggles[6*i+1] = er2;
        wiggles[6*i+2] = er3;

        wiggles[6*i+3] = et1;
        wiggles[6*i+4] = et2;
        wiggles[6*i+5] = et3;
    }

    //console.log(max_t, max_ti, max_r, max_ri)
}

function get_timestamp() {
    return (new Date()).toISOString().substr(0, 19);
}

function get_axis(r1, r2, r3) {
    rnorm = Math.sqrt(r1 * r1 + r2 * r2 + r3 * r3)
    k1 = r1 / rnorm
    k2 = r2 / rnorm
    k3 = r3 / rnorm
    return new THREE.Vector3(k1, k2, k3)
}

function get_angle(r1, r2, r3) {
    return Math.sqrt(r1 * r1 + r2 * r2 + r3 * r3)
}

/* get position of a camera center, which is -R'*t */
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
        CURRENT_EIG -= 1
        if ( CURRENT_EIG < 0 ) {
            CURRENT_EIG = n_eigs - 1;
        }
        set_eigenvector(CURRENT_EIG);
    } else if ( e.key == 'R' ) {
        //REC = true;
        capturer.start();
    } else if ( e.key == 'S' ) {
        //REC = false;
        //capturer.stop();
        //capturer.save();
        capturer = new CCapture( {
            format: 'png',
            framerate: 120,
            verbose: false,
            timeLimit: 1,
            autoSaveTime: 10,
            name: get_timestamp()
        } );
    } else if ( e.key == 'V' ) {
        // toggle static mode (static frusta, camera center points)
        // vs dynamic mode (animated frusta)
        STATIC = !STATIC 
        camera_points.visible = !camera_points.visible
    } else if ( e.key == 'p' ) {
        scene_points.visible = !scene_points.visible // toggle scene point visibility
    } else if ( e.key == "+" ) {
        CAM_PT_SCALE *= 1.1;
    } else if ( e.key == "_" ) {
        CAM_PT_SCALE *= 1/1.1;
    } else if ( e.key == '?' ) {
        var overlay = document.getElementById("kb_overlay")
        if (overlay.style.visibility == "hidden") {
            overlay.style.visibility = "visible"
        } else {
            overlay.style.visibility = "hidden"
        }
    }
    update_overlay()
}

function update_overlay() {
    var overlay = document.getElementById("overlay")
    display_html = "<table>";
    display_html += table_row("Eigenvector:", CURRENT_EIG)
    display_html += table_row("Eigenvalue:", CURRENT_EV)
    //display_html += table_row("Recording:", REC)
    display_html += table_row("Eig T scale:", EIG_T_FUDGE)
    display_html += table_row("Eig R scale:", EIG_R_FUDGE)
    display_html += "</table>"
    overlay.innerHTML = display_html
}

function table_row(label, value) {
    return "<tr><td>" + label + "</td><td>" + value + "</td></tr>";
}

// Update function to wiggle points by their wiggle vectors * sin(t)
// this is seemingly CPU-based but the framerate only seems to slow
// down somewhere between 100,000 and 1 million particles.
function wiggle_cameras(t) {
    var buffer_pos = camera_points.geometry.attributes.position.array;

    if (!STATIC) {
        sint = Math.sin(t);
    } else {
        sint = Math.abs((t % 4)-2) - 1

        // vary the color of the point to show phase
        color_t = Math.abs((t % 4) - 2) / 2

        r = Math.max(color_t, 0.5)
        g = Math.max(-1 * color_t + 0.5, 0)
        b = Math.max(1-color_t, 0.5)
        camera_points.material.color.setRGB(r, 0.5, b)

        // make point shrink toward the edge for slightly pointy lines
        camera_points.material.size = CAM_PT_SCALE * (-Math.abs((t*2) % 4 - 2)/4 + 1)
    }

    for ( var i = 0; i < n_cameras; i++ ) {

        var cam = frusta[i];

        var x = positions[3*i] + EIG_T_FUDGE * wiggles[6*i+3] * sint;
        var y = positions[3*i+1] + EIG_T_FUDGE * wiggles[6*i+4] * sint;
        var z = positions[3*i+2] + EIG_T_FUDGE * wiggles[6*i+5] * sint;

        if ( STATIC ) {
            buffer_pos[3*i] = positions[3*i] + EIG_T_FUDGE * wiggles[6*i+3] * sint;
            buffer_pos[3*i+1] = positions[3*i+1] + EIG_T_FUDGE * wiggles[6*i+4] * sint;
            buffer_pos[3*i+2] = positions[3*i+2] + EIG_T_FUDGE * wiggles[6*i+5] * sint;

            cam.position.set(positions[3*i], positions[3*i+1], positions[3*i+2])

        } else {
            var r1 = rotations[3*i] + EIG_R_FUDGE * wiggles[6*i] * sint;
            var r2 = rotations[3*i+1] + EIG_R_FUDGE * wiggles[6*i+1] * sint;
            var r3 = rotations[3*i+2] + EIG_R_FUDGE * wiggles[6*i+2] * sint;

            axis = get_axis(r1, r2, r3);
            angle = get_angle(r1, r2, r3);

            cam.position.set(x, y, z);
            cam.setRotationFromAxisAngle(axis, angle)
        }


    }

    // need to do this after changing positions
    //camera_points.geometry.computeBoundingSphere();
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

    capturer.capture( renderer.domElement )
    wiggle_cameras( t_elapsed * 4 );
    geometry.attributes.position.needsUpdate = true;

};

render();
