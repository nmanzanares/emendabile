//import "geometry.js";

let width = 0;
let height = 0;
let qvga = {width: {exact: 320}, height: {exact: 240}};
let vga = {width: {exact: 640}, height: {exact: 480}};
let resolution = window.innerWidth < 640 ? qvga : vga;

let param_adaptive_threshold_block_size = 45;
let param_adaptive_threshold_offset = 4;

let param_hough_thresholds = [180, 160, 140, 120, 100, 80];
let hough_threshold = 30;

let th_filter_rho =6;
let param_directions_threshold = 0.4;
let param_check_corners_tolerance_mul = 6;

//DECISION
let decisions;
let puntation;
let answer_cells;
let corners;
let type;
let param_conter=0;

param_cross_mask_margin = 0.6;
param_cross_mask_margin_2 = 0.75;
param_cross_mask_thickness = 0.2;
param_cross_mask_threshold = 0.08;

function process(src){
  changeCam();
  let pre = pre_process(src);
  console.log("param: " + param_hough_thresholds[param_conter]);
  let lines_ob_mat = detect_lines(pre, param_conter);
  let nf_lines = draw_lines(src, lines_ob_mat);
  console.log(lines_ob_mat.rows);
  if(lines_ob_mat.rows==0){
  	return false;
  }
  if(debug_mode){
  	cv.imshow('canvasOutput', nf_lines);
  }
  
  let lines_nf = mat_to_array(lines_ob_mat);
  lines_nf.sort(function(a, b){return a[1] - b[1]}); //Ordenado por theta
  
  let axes = detect_boxes(lines_nf,model.dims);
  
  
  if(axes!=null && axes.length==2){
  	axes = filter_axes(axes, 6);
  	axes = filter_axes(axes, 12);
  	let f_lines = draw_axes(src, axes);
	
	//cv.imshow('canvasOutput', f_lines);
	
	corners = cell_corners(axes[1].lines, axes[0].lines, src.cols, src.rows, model.dims);
	if(corners.length == 0){
		console.log('Esquinas no detectadas');
		return false;
	}	
	let lines_points = draw_corners(f_lines,corners);
	
	decisions = decide_cells(corners, pre);
	type = decide_type(corners, pre);
	console.log("Tipo de examen: " + type);
	if(type==-1){
		return false;
	}
	if(type>model.sols.length-1){
		final_div.insertAdjacentHTML('afterbegin', '<span class="error"> Por favor, introduzca enel modelo soluciones para el tipo ' + (type+1) + ' </span>');
	}
	draw_decisions(src,model.sols,type,decisions,corners);
 
	//cv.imshow('canvasOutput', lines_points);
  }else{
  	console.log('Se han obtenido unos ejes nulos');
  	return false;
  }
  
  puntation = correction(model.sols,type,decisions,model.plus,model.deduct,model.emptiness);
  print_puntation(puntation);
  
  pre.delete();
  nf_lines.delete();
  lines_ob_mat.delete();
  
  return true;
}

function findIdVideo(){
    return new Promise(function(resolve, reject) {
	navigator.mediaDevices.enumerateDevices()
	    .then(function(devices) {
                var ids = new Array();
		devices.forEach(function(device) {
		    if(device.kind === "videoinput"){
                        ids.push(device.deviceId);
                    }
		});
                resolve(ids);
	});
    });
}

let video = document.getElementById("video");
let stream = null;
let vc = null;
let streaming = false;
let deviceIds;
const FPS = 50;
let cameras = new Array();
let currentCamera;

findIdVideo().then(function(cams) {
    cameras = cams;
    currentCamera = 0;
});

function startCamera(){
    startAndStop.removeAttribute('hidden');
    finish_button.removeAttribute('hidden');
    startAndStop.innerText = 'Stop Video';
    canvasOutput.setAttribute('hidden','');
    video.removeAttribute('hidden');
    if(streaming) return;
    navigator.mediaDevices.getUserMedia(
        {video: {
            deviceId: {
                exact: cameras[currentCamera]
            }
        },
         audio: false})
        .then(function(s) {
            stream = s;
            video.srcObject = s;
            video.play();
        })
    .catch(function(err) {
    console.log("An error occured! " + err);
  });

  video.addEventListener("canplay", onVideoCanPlay, false);
}

function changeCam(){
	if(cameras.length > 1 && streaming){
		changeCamera.removeAttribute('hidden');
	}else{
		changeCamera.setAttribute('hidden','');
	}
	if(currentCamera> cameras.length-1){
		currentCamera=0;
	}
}


function startVideoProcessing(){
	if (!streaming) { console.warn("Please startup your webcam"); return; }
	src = new cv.Mat(height, width, cv.CV_8UC4);
	let begin = Date.now();
 	vc.read(src);
 	success = process(src);
	if(!success){
		let delay = 1000/FPS - (Date.now() - begin);
        setTimeout(startVideoProcessing, delay);
        param_conter++;
        if(param_conter== param_hough_thresholds.length) param_conter=0;
	}else{
		stopCamera();
	}
}

function stopCamera(){
	if (this.video) {
       video.pause();
       video.srcObject = null;
       video.removeEventListener('canplay', onVideoCanPlay);
    }
    if (stream) {
       stream.getVideoTracks()[0].stop();
    }
    streaming=false;

	video.setAttribute('hidden','');
	canvasOutput.removeAttribute('hidden');
    startAndStop.innerText = 'Next Exam';
}

function onVideoCanPlay (ev){
    if (!streaming) {
	  width=600;
	  height=400;
      video.setAttribute("width", width);
      video.setAttribute("height", height);
      streaming = true;
      vc = new cv.VideoCapture(video);
    }
    startVideoProcessing();
  }

function pre_process(src) { //1
  let dst = new cv.Mat();
  cv.cvtColor(src, dst, cv.COLOR_RGBA2GRAY);
  cv.adaptiveThreshold(dst, dst, 200, cv.ADAPTIVE_THRESH_GAUSSIAN_C, cv.THRESH_BINARY_INV, param_adaptive_threshold_block_size, param_adaptive_threshold_offset);
  //cv.Canny(dst, dst, 50, 100, 3, false);
  return dst;
}

function canny(src){ //1
	let dst = new cv.Mat();
	cv.cvtColor(src, src, cv.COLOR_RGBA2GRAY, 0);
	cv.Canny(src, dst, 50, 100, 3, false);
	return dst;
}

function detect_lines(image, param){ //2
	let lines_mat = new cv.Mat();
	//let j=0;
	//let mos=0;
	//let mos2=0;
	/*do{
		cv.HoughLines(image, lines_mat, 1, Math.PI / 180,
		          param_hough_thresholds[j], 0, 0, 0, Math.PI);
		if(lines_mat.rows > mos && lines_mat.rows < 100){
			mos = lines_mat.rows;
			mos2 = j;
		}
		j++;
	} while (j<param_hough_thresholds.length)
	*/
	cv.HoughLines(image, lines_mat, 1, Math.PI / 180,
		          param_hough_thresholds[param], 0, 0, 0, Math.PI);
		          	         
	return lines_mat;	
}

function mat_to_array(mat){ //3
	arr = new Array();
	for(let i = 0; i < mat.rows; ++i) {
		el = [mat.data32F[i*2], mat.data32F[i*2 + 1]];
		arr.push(el);
	}
	return arr;
}

function expected_lines(dimensions){ 
	let expected = {h:0, v:0};
	let c=0;
	let c_max=0;
	for(let i=0;i<dimensions.length;i++){c+=dimensions[i][0];} 
	expected.v = dimensions.length + c;
	for(let i=0;i<dimensions.length;i++){
		c=dimensions[i][1];
		if(c>c_max){c_max=c;}
	} 
    expected.h = 1 + c_max;
    return expected;
}

function detect_directions(lines){ //5
	let axes = new Array();
	let rho = lines[0][0];
	let theta = lines[0][1];
	axes.push({theta:theta, lines:[[rho, theta]]});
	for(let i=1; i<lines.length;i++){
		rho = lines[i][0];
		theta = lines[i][1];
	    if(Math.abs(theta - axes[axes.length-1].theta) < param_directions_threshold){
	        axes[axes.length-1].lines.push([rho, theta])
	    }else{
			axes.push({theta:theta, lines:[[rho, theta]]});
	    }
	}
	//Juntamos las líneas en ejes separados a distancia pi
	if(Math.abs(axes[0].theta - axes[axes.length-1].theta + 3.1416) < param_directions_threshold){
		let ax = axes.pop();
		for(let i=0; i<ax.lines.length;i++){
			axes[0].lines.push([-ax.lines[i][0], ax.lines[i][1] - 3.1416]);
		}
	} 
	//La theta de referencia del eje será la media de todas sus thetas       
	for(let i=0; i < axes.length; i++){
		let sum_theta=0;
		for(let j=0; j<axes[i].lines.length; j++){
			sum_theta += axes[i].lines[j][1];
		}
		avg = sum_theta/axes[i].lines.length;
	    axes[i].theta = avg;
		axes[i].lines.sort(function(a, b){return a[0] - b[0]})
	}
	//Colocamos en primer lugar el eje vertical
	if(Math.abs(axes[axes.length-1].theta - 3.1416) < Math.abs(axes[0].theta)){
	    let aux = axes.splice(0,1,axes[axes.length-1]);
	    axes.pop();
	    axes.push(aux.pop());
	}
	return axes
}

function detect_boxes(lines, dimensions){ //5
	let expected = expected_lines(dimensions);
	let axes = detect_directions(lines);
	for(let i=0; i<axes.length; i++){
		if(axes[i].lines.length < Math.min(expected.h, expected.v)){
			axes.splice(i,1);
			i--;
		}
	}
	if(axes.length == 3){
		if(angles_perpendicular(axes[0].theta,axes[1].theta)){
			axes.splice(2,1);
		}else{
			if(angles_perpendicular(axes[0].theta,axes[2].theta)){
				axes.splice(1,1);
			}else{
				axes.splice(0,1);
			}
		}
	}else{
		if(axes.length == 4){
			for(let j=0; j<axes.length; j++){
				if(Math.abs(axes[j].theta)>param_directions_threshold && Math.abs(axes[j].theta)< (3.1416 - param_directions_threshold) && Math.abs(axes[j].theta - (3.1416/2)) > param_directions_threshold){
					axes.splice(j,1);
					j--;
				}
			}
		}
	}
	if(axes.length == 2){
		if(angles_perpendicular(axes[0].theta,axes[1].theta)){
			return axes;
		}else{
			return null;
		}
	}else{
		return null;
	}
}

function filter_axes(axes, param){ //6
	let hlines = collapse_lines(axes[1].lines, param);
	let vlines = collapse_lines(axes[0].lines, param);
	if(hlines==null || vlines==null){
		return axes;
	}else{
		let axis = new Array();
		axis.push({theta:axes[0].theta, lines: vlines});
		axis.push({theta:axes[1].theta, lines: hlines});
		return axis;
	}
}

function collapse_lines(lines_nf, param){ //6
	lines_nf.sort(function(a, b){return a[0] - b[0]}) //Ordenado por rho
 	let sum_rho = lines_nf[0][0];
	let sum_theta = lines_nf[0][1];
    let num_lines = 1;
    let last_line = lines_nf[0];
	let lines = new Array();
	
	for (let i = 1; i < lines_nf.length; ++i){
		let rho = lines_nf[i][0];
		let theta = lines_nf[i][1];		
		if((Math.abs(rho-last_line[0]) > param)){
			lines.push([sum_rho/num_lines,sum_theta/num_lines]);
			sum_rho = rho;
			sum_theta = theta;
			num_lines = 1;
			last_line = lines_nf[i];
		}else{
			sum_rho += rho;
			sum_theta += theta;
			num_lines++;
			last_line = lines_nf[i];
		}
	}
	lines.push([sum_rho/num_lines,sum_theta/num_lines]);
	return lines;	
}

function check_corners(corner_matrixes, width, height){
	//Comprobamos distancias entre las líneas horizontales
    let corners = corner_matrixes[0];
    let ycorn = new Array();
    for(let i=0; i < corners.length; i++){
    	let row = corners[i];
    	ycorn.push(row[row.length-1].y);
    }
    let dist = new Array();
    let diff_dist = new Array();
    for(let i=1; i<ycorn.length; i++){
    	dist.push(ycorn[i]-ycorn[i-1]);
    }
    for(let i=1; i<dist.length; i++){
    	diff_dist.push(dist[i]-dist[i-1]);
    }
    let max_diff = 1 + ((Math.max.apply(null,dist) - Math.min.apply(null,dist)) / dist.length) * param_check_corners_tolerance_mul;
	if(Math.max.apply(null,diff_dist) > max_diff){
		return false;
	}
	if(0.5*Math.max.apply(null,dist) > Math.min.apply(null,dist)){
		return false;
	}

	for(let i=0; i<corner_matrixes.length; i++){
		let corners = corner_matrixes[i];
		for(let j=0; j<corners.length; j++){
			let row = corners[j];
			for(let k=0; k< row.length; k++){
				let p = row[k];
				if(p.x <0 || p.x > width || p.y<0 || p.y>height){
					return false;
				}
			}
		}
	}
	
	for(let i=0; i<corner_matrixes.length; i++){
		let corners = corner_matrixes[i];
		for(let j=0; j<corners.length-1; j++){
			for(let k=0; k< corners[0].length-1; k++){
				if(corners[j][k].y >= corners[j+1][k].y || corners[j][k+1].y >= corners[j+1][k+1].y || corners[j][k].x >= corners[j][k+1].x || corners[j+1][k].x >= corners[j+1][k+1].x){
					return false;
				}
			}
		}
	}
	
	return true;
}

function cell_corners(hlines, vlines, src_w, src_h, dimensions){ //7
	let expected = expected_lines(dimensions);
	let corner_matrixes = new Array();
	if(vlines.length != expected.v){
        if(vlines.length > expected.v && vlines.length <= expected.v + 2){
            // Remove one or two spurious lines
            // vlines = discard_spurious_lines(vlines, v_expected)
	 		return corner_matrixes;
        }else{
        	return corner_matrixes;
        }
    }
    if(hlines.length < expected.h){
    	return corner_matrixes;
    }else{
    	if(hlines.length > expected.h){
			hlines1 = hlines.slice(hlines.length - expected.h, hlines.length);
			hlines2 = hlines.slice(0,hlines.length - expected.h);
			hlines = hlines1.concat(hlines2);		
    	}
    }
    
    let v_i= 0;
   	for(let i=0; i<dimensions.length; i++){
   		let w= dimensions[i][0];
		let h= dimensions[i][1];
		let corners = new Array();
    	for(let j=0; j<h+1; j++){
    		let p_group = new Array();
    		for(let k=v_i; k<v_i+w+1; k++){
    			let p = intersection(hlines[j],vlines[k]);
    			p_group.push(p);
    		}
    		corners.push(p_group);
    	}
    	corner_matrixes.push(corners);
    	v_i += 1+w;
    }
    
    if(check_corners(corner_matrixes,src_w, src_h)){
    	return corner_matrixes;
    }else{
    	return new Array();
    }	
}

function answer_cells_geometry(corners){
	let cells = new Array();
	for(let i = 0; i < corners.length; ++i){
		for(let j = 0; j < corners[i].length-1; j++) {
			let row = new Array();
			for(let k = 0; k < corners[i][j].length-1; k++){
				let cell = {plu: corners[i][j][k], pru: corners[i][j][k+1],pld: corners[i][j+1][k], prd: corners[i][j+1][k+1]};
				row.push(cell);
			}
			cells.push(row);
		}
	}
	return cells;
}

function is_cross(pre, mask, masked, plu, pru, pld, prd){
	let thickness = distance(plu,pru) * param_cross_mask_thickness;
	let closer_corners_1 = closer_points_rel(plu, prd, param_cross_mask_margin,thickness/2);
	let closer_corners_2 = closer_points_rel(pru, pld, param_cross_mask_margin,thickness/2);
	draw_cross_mask(mask, closer_corners_1, closer_corners_2, thickness);
	
	closer_corners_1 = closer_points_rel(plu, prd, param_cross_mask_margin_2,thickness/4);
	closer_corners_2 = closer_points_rel(pru, pld, param_cross_mask_margin_2,thickness/4);
	draw_cross_mask(mask, closer_corners_1, closer_corners_2, thickness/2);

	let mask_pixels = cv.countNonZero(mask);
	cv.multiply(pre,mask,masked);
	let masked_pixels = cv.countNonZero(masked);
	cv.multiply(cv.Mat.zeros(mask.size().height, mask.size().width, mask.type()),mask,mask);
	cv.multiply(cv.Mat.zeros(masked.size().height, masked.size().width, masked.type()),masked,masked);
	let cell_marked;
	if(masked_pixels > param_cross_mask_threshold * mask_pixels){
		cell_marked=1;
	}else{
		cell_marked=0;
	}
	return cell_marked;
}

function decide_cells(corner_matrixes, pre){
	answer_cells = answer_cells_geometry(corner_matrixes);
	let decisions = new Array();
	let mask = new cv.Mat(pre.size().height, pre.size().width, pre.type());
	let masked = new cv.Mat(pre.size().height, pre.size().width, pre.type());
	
	for(let i = 0; i < answer_cells.length; ++i){
		let row_decisions= new Array();
		for(let j = 0; j < answer_cells[i].length; j++) {
			let cell = answer_cells[i][j];
			let decision = is_cross(pre, mask, masked, cell.plu, cell.pru, cell.pld, cell.prd);
			row_decisions.push(decision);		
		}
		decisions.push(row_decisions);
	}
	return decisions;
}

function decide_type(corners, pre){
	answer_cells = answer_cells_geometry(corners);
	let bits = new Array();
	let mask = new cv.Mat(pre.size().height, pre.size().width, pre.type());
	let masked = new cv.Mat(pre.size().height, pre.size().width, pre.type());
	
	last1 = Math.round(answer_cells.length/2)-1;
	last2 = Math.trunc(answer_cells.length/2)-1;
	
	for(let i=0; i<answer_cells[last1].length; i++){
		let centers = get_centers(answer_cells[last1][i]);
		bits.push(decide_bits(pre, centers, mask, masked));
			
	}
	
	let error_ib;
	let type = bits[0] + bits[1]*2 + bits[2]*4;
	for(let j=0; j<4; j++){
		if(j+4 < bits.length || bits[j]==-1){
			if(bits[j] != bits[j+4] || bits[j]==-1){
				error_ib = true;
			}
		}
	}
	if(error_ib){
		console.log("Error en la lectura de infobits");
		return -1;
	}else{
		return type;
	}
}

function decide_bits(pre, centers, mask, masked){
	cv.multiply(cv.Mat.zeros(mask.size().height, mask.size().width, mask.type()),mask,mask);
	cv.multiply(cv.Mat.zeros(masked.size().height, masked.size().width, masked.type()),masked,masked);
	
	mask_d = mask.clone();
	masked_d = masked.clone();
	
	let color = new cv.Scalar(255, 0, 0);
	let center_up = new cv.Point(centers.up.x, centers.up.y);
	cv.circle(mask, center_up, centers.r, color, -1);
	let mask_pixels_up = cv.countNonZero(mask);
	cv.multiply(pre,mask,masked);
	let masked_pixels_up = cv.countNonZero(masked);
	
	let center_down = new cv.Point(centers.down.x, centers.down.y);
	cv.circle(mask_d, center_down, centers.r, color, -1);
	let mask_pixels_down = cv.countNonZero(mask_d);
	cv.multiply(pre,mask_d,masked_d);
	let masked_pixels_down = cv.countNonZero(masked_d);
	
	let bit;
	if(masked_pixels_up > param_cross_mask_threshold * mask_pixels_up){
		bit=1;
	}else{
		if(masked_pixels_down > param_cross_mask_threshold * mask_pixels_down){
			bit=0;
		}else{
			bit=-1;
		}
	}
	return bit;
}

function get_centers(cell){
	let center_up = {x: (cell.pld.x + cell.prd.x)/2, y: (cell.pld.y + (cell.pld.y - cell.plu.y)/2 )};
	let center_down = {x: (cell.pld.x + cell.prd.x)/2, y: (cell.pld.y + (cell.pld.y - cell.plu.y)*3/2 )};

	return {up: center_up, down: center_down, r: Math.abs((cell.pld.x - cell.prd.x)/4)};
}

function is_clicked(cX,cY,cell){
	if(cX > cell.plu.x && cX < cell.pru.x && cY > cell.plu.y && cY < cell.pld.y){
		return true;
	}else{
		return false;
	}
}

function correction(sols, type, decisions, plus, deduct, emptiness){
	let correct=0;
	let incorrect=0;
	let empty=0;
	let row_marked=false;
	let puntuation;
	let solutions = sols[type];
	for (let i = 0; i < decisions.length; i++) {
		for (let j = 0; j < decisions[i].length; j++) {
			if(decisions[i][j]==1){
				row_marked=true;
				if(solutions[i][j]==1){
					correct++;
				}else{
					incorrect++;
				}	
			}
		}
		if(!row_marked){
			empty++;
		}
		row_marked=false;
	}
	puntuation = {note:(plus*correct + deduct*incorrect + emptiness*empty), correct: correct, incorrect:incorrect, empty:empty};
	
	return puntuation; 
}

function check_click(src,cX,cY){
	console.log("cX: " + cX + ", cY: " + cY);	
	for(let i = 0; i < answer_cells.length; ++i){
		for(let j = 0; j < answer_cells[i].length; j++) {
			let cell = answer_cells[i][j];
			if(is_clicked(cX,cY,cell)){
				console.log("Celda: plu: " + cell.plu.x + "," + cell.plu.y + " pru: " + cell.pru.x + "," + cell.pru.y + " pld: " + cell.pld.x + "," + cell.pld.y + " prd: " + cell.prd.x + "," + cell.prd.y);
				if(decisions[i][j]==1){
					decisions[i][j]=0;
					console.log("Ponemos la celda " + (i+1) + "," + (j+1) + " vacía");
				}else{
					for(let k = 0; k < answer_cells[i].length; k++) {
						if(k==j){
							decisions[i][k] = 1;
							console.log("Ponemos la celda " + (i+1) + "," + (k+1) + " marcada");
						}else{
							decisions[i][k] = 0;
						}
					}
				}
			}
		}
	}
	draw_decisions(src,model.sols,type,decisions,corners);
	puntation = correction(model.sols, type, decisions, model.plus, model.deduct, model.emptiness);
	print_puntation(puntation);
}

function draw_cross_mask(mask, closer_corners_1, closer_corners_2, thickness){ //OJO
	cv.line(mask, closer_corners_1[0], closer_corners_1[1], [255,255,255,0], Math.round(thickness));
	cv.line(mask, closer_corners_2[0], closer_corners_2[1], [255,255,255,0], Math.round(thickness));
}

function draw_decisions(src, sols, type, decisions, corners){ //10
	let dst = src.clone();
	let solutions = sols[type];
	let answer_cells = answer_cells_geometry(corners);
	for (let i = 0; i < answer_cells.length; i++) {
		for (let j = 0; j < answer_cells[i].length; j++) {
			if(decisions[i][j]==1){
				let x0 = mean_corn(answer_cells[i][j].plu.x, answer_cells[i][j].pru.x, answer_cells[i][j].pld.x, answer_cells[i][j].prd.x);
				let y0 = mean_corn(answer_cells[i][j].plu.y, answer_cells[i][j].pru.y, answer_cells[i][j].pld.y, answer_cells[i][j].prd.y);
				let center = {x:x0, y:y0};
				let radio = (answer_cells[i][j].pru.x - answer_cells[i][j].plu.x)/2
				if(solutions[i][j]==1){
					cv.circle(dst, center, radio, [0,255,0,255]);
				}else{
					cv.circle(dst, center, radio, [255,0,0,255]);
				}	
			}
		}
	}
	cv.imshow('canvasOutput', dst);
}

function draw_lines(src, lines){ //4
	let dst = src.clone();
	for (let i = 0; i < lines.rows; ++i) {
		let rho = lines.data32F[i * 2];
		let theta = lines.data32F[i * 2 + 1];
		let a = Math.cos(theta);
		let b = Math.sin(theta);
		let x0 = a * rho;
		let y0 = b * rho;
		let startPoint = {x: x0 - 1000 * b, y: y0 + 1000 * a};
		let endPoint = {x: x0 + 1000 * b, y: y0 - 1000 * a};
		cv.line(dst, startPoint, endPoint, [255, 0, 0, 255]);
	}
	return dst;
}

function draw_corners(src, corner_matrixes){ //8
	let dst = src.clone();
	for (let i = 0; i < corner_matrixes.length; i++) {
		for (let j = 0; j < corner_matrixes[i].length; j++) {
			for (let k = 0; k < corner_matrixes[i][j].length; k++) {
				let x0 = corner_matrixes[i][j][k].x;
				let y0 = corner_matrixes[i][j][k].y;
				let center = {x:x0, y:y0};
				cv.circle(dst, center, 5, [0, 0, 255, 255]);
			}
		}
	}
	return dst;
}

function draw_lines_arr(src, lines){ //4
	let dst = src.clone();
	for (let i = 0; i < lines.length; ++i) {
		let rho = lines[i][0];
		let theta = lines[i][1];
		let a = Math.cos(theta);
		let b = Math.sin(theta);
		let x0 = a * rho;
		let y0 = b * rho;
		let startPoint = {x: x0 - 1000 * b, y: y0 + 1000 * a};
		let endPoint = {x: x0 + 1000 * b, y: y0 - 1000 * a};
		cv.line(dst, startPoint, endPoint, [255, 0, 0, 255]);
	}
	return dst;
}

function draw_axes(src, axes){ //4
	let dst = src.clone();
	for(let i=0; i<axes.length;i++){
		for (let j = 0; j < axes[i].lines.length; j++) {
			let rho = axes[i].lines[j][0];
			let theta = axes[i].lines[j][1];
			let a = Math.cos(theta);
			let b = Math.sin(theta);
			let x0 = a * rho;
			let y0 = b * rho;
			let startPoint = {x: x0 - 1000 * b, y: y0 + 1000 * a};
			let endPoint = {x: x0 + 1000 * b, y: y0 - 1000 * a};
			cv.line(dst, startPoint, endPoint, [255, 0, 0, 255]);
		}
	}
	return dst;
}





