var cell = 10;//mm
var L = 100;//mm
//var height = 40;//cells
//var width = 40;//cells
var v = 4000; //mm/s

var c_metal = 237 /1000; //aluminum
var c_air = 0.25 /1000;
var c_contact = 50 / 1000/1000; //guesswork, ~5-10 for steel, rough surface = more?

var max = Math.max;
var min = Math.min;

var GlobalState={'flow':[],'heat':[]}

/*\
      frame := a 4x1 array;
    starting in X-direction (1,0)
    containing values relevant to
    the EDGES of evaluation-window
\*/

/* note: all following constants have a more complete form, noted below.
   however, since internal (metal to metal, air to air) is normed to 1
   we are only interested in relative conduction differences.
   
   metal_c := cell*L / (cell/2) * 237
   air_c := cell*L / (cell/2) * 0.25
   cell*L is the conduction cross section
   cell/2 is the distance between cell-center and wall
   237 is conduction in aluminum
   0.25 is conduction in air (NOT convection!)

   gap_c = 7 * cell*L
   7 is the conduction approximation for metal-to-air interface
   (therefore only area, not depth, is needed for this calculation)

   total conduction across barrier is given by:
        metal_c*air_c*gap_c/(metal_c+air_c+gap_c)
   will be normalised to (=divided by) metal_c and air_c, e.g.:*/
// metal_to_air = air_c*gap_c/(metal_c+air_c+gap_c)
   var metal_to_air = c_air * c_contact * cell*L /(c_metal + c_air + c_contact * cell/2)
// and similarly:
   var air_to_metal = c_metal * c_contact * cell*L /(c_metal + c_air + c_contact * cell/2)

/* when calculating an air-cell, coupling with fresh air has to be considered.
   This is achieved by considering the cell as voxl in contact with 2
   "virtual" cells at ambient temperature. The heat conducted
   can be calculated as the heat transported by the moving volume
   thanks to it's heat capacity. Conductance therefore depends on velocity
   which is stored in the flow-matrix for each cell. Normalisation is done
   by dividing with air-air-conductance as it isn't relevant in metal-cells.
   Since hot air is leaving the cell and cold air is entering, the cell is
   indeed in contact with 2 cells. This is accomplished by multiplying with
   a factor of 2.*/
   var extra_air = v * 2 * cell**2 * 1.01 * 1.15 / (L*c_air);
//             v * 2 *    A    * c_p  *  rho / air_air_c

console.log(air_to_metal);
console.log(metal_to_air);
console.log(extra_air);
function get_fresh_array(w,h,init=0){
	var ret = new Array(w);
	for(let i=0; i<w; i++){
		var colm = new Float32Array(h);
		//we don't actually need initialisation ...
		//for(let j=0; j<h; j++) colm[j] = init;
		ret[i] = colm;
	}
	return ret;
}

function fresh_from_copy(orig, init=0){
	return get_fresh_array(orig.length, orig[0].length, init);
}

function prepare_air(width, height, source){
	//carry around 1 unnecessarry row
	var air  = get_fresh_array(width+2, height+2);
	/*for (var j=0; j<height+2; j++){
		air[0][j]=1
		air[1][j]=1
		air[width][j]=1
		air[width+1][j]=1
	}*/
	for(var i=0; i<air.length; i++){
		//highest row mirrors/isolates border condition
		//air[i][height]=height*air_stepsize
		//highest row is moving "plane"
		air[i][height+1] = 1;
		//lowest row is ref, j=1 is full of metal
		for (var j=0; j<height+2; j++){
			air[i][j] = (j+1)/(height+1);
		}
		for(let i=source[0];i<(source[1]+1);i++){
			air[i][0]=0;
			air[i][1]=0;
		}
		//air[i][0]= (i>2&&i<air.length-3) ? 0 : 0.1
		//air[i][1]= (i>2&&i<air.length-3) ? 0 : 0.1
	}
	return air
}

function prepare_fire(width, height, source){
	var fire = get_fresh_array(width+2, height+2);
	for (var i=1; i<fire.length-2; i++){
		if(i<source[0]+1){
			fire[i][1]=(1/2)//*i/source[0];
			//isolation border condition
			fire[i][0]=(1/2)//*i/source[0];
		} else if(i>source[1]){
			fire[i][1]=(1/2)//*(width-i+1)/(width-source[1]+1);
			//isolation border condition
			fire[i][0]=(1/2)//*(width-i+1)/(width-source[1]+1);
		} else {
			fire[i][1]=0.5;
			//force heat source high
			fire[i][0]=1;
		}
	}
	//isolation border conditions
	      fire[0][1]=(1/2)/source[0];
	      fire[0][0]=(1/2)/source[0];
	fire[width+1][1]=(1/2)/(width-source[1]);
	fire[width+1][0]=(1/2)/(width-source[1]);
	//technically, corner-cells won't ever be touched, but ... futur-proof or something.
	return fire;
}


/*\
| | iterstep
\*/
function step_flow(flow, error_as_array=false){
	var new_flow=fresh_from_copy(flow);
	var err = error_as_array? []:0;
	//don't touch frame (left-/right-most column, top/bottom row)
	for (let i=1; i<flow.length-1; i++){
		for (let j=1; j<flow[0].length-1; j++){
			//don't change if already heatsink!
			if(flow[i][j] == 0) continue;

			new_flow[i][j] = (flow[i-1][j]+flow[i][j-1]+flow[i][j+1]+flow[i+1][j])/4.0;
			if (error_as_array){
				err.append((new_flow[i][j] - flow[i][j])**2)
			} else {
				err +=  (new_flow[i][j] - flow[i][j])**2
			}
		}
	}
	return {"result":new_flow,"error":err}
}

function step_heat(heat, flow, error_as_array=false){
	var new_heat=fresh_from_copy(flow);
	var err= error_as_array? [] : 0;
	//don't touch outer lines, handle later
	for (let i=1; i<heat.length-1; i++){
		//don't touch upper most row, handle later
		for (let j=1; j<heat[0].length-1; j++){
			//amiheatsink?
			let amsink= (flow[i][j] == 0);

			//adding extra_air to weight. 
			var weight = amsink? 0 : extra_air*flow[i][j];
			//no extra_air for sum, as it's value is 0 (=ambient temp)
			var sum = 0;
			//circular neighborhood
			let neighbors = [
				{'i':i,'j':j-1},
				{'i':i-1,'j':j},
				{'i':i+1,'j':j},
				{'i':i,'j':j+1}
			];
			for (let n=0;n<4;n++){
				//shorthand, o bc. circular
				let o=neighbors[n];
				let issink = (flow[o.i][o.j] == 0);

				if ((amsink && issink)||(!amsink && !issink)){
					sum += heat[o.i][o.j];
					weight += 1;
				} else if (amsink && !issink){
					sum += heat[o.i][o.j] * metal_to_air;
					weight += metal_to_air;
				} else if (!amsink && issink){
					sum += heat[o.i][o.j] * air_to_metal;
					weight += air_to_metal;
				}
			}
			new_heat[i][j] = sum? sum/weight : 0;
			if (error_as_array){
				err.push((new_heat[i][j] - heat[i][j])**2);
			} else {
				err +=  (new_heat[i][j] - heat[i][j])**2;
			}
		}
	}
	return {'result':new_heat, 'error':err};
}

//!!! needs equilibrated flow !!!
function equilibrate(grid, neumann, flow=false, max_error=0.1, max_loop=100){
	var new_grid={'result':grid};
	for(let totally_unnecessarry_counter=0; totally_unnecessarry_counter < max_loop; totally_unnecessarry_counter++){
		if(flow){
			new_grid=step_heat(new_grid.result, flow);
		} else {
			new_grid=step_flow(new_grid.result)
		}
		//apply neumann (left and right) by copying inner columns
		if (neumann[0]){
			new_grid.result[0]             = new_grid.result[1];
		} else {
			new_grid.result[0]             = grid[0];
		}
		if (neumann[2]){
			new_grid.result[grid.length-1] = new_grid.result[grid.length-2];
		} else {
			new_grid.result[grid.length-1] = grid[grid.length-1];
		}
		//apply neumann/dirichlet (top and bottom) by looping
		for (let i=0; i<grid.length; i++){
			if (neumann[1]){
				//if source
				if (grid[i][grid[0].length-1] == 1){
					//apply dirichlet
					new_grid.result[i][grid[0].length-1] = 1
				//else (=not source)
				} else {
					//apply neumann
					new_grid.result[i][grid[0].length-1] = new_grid.result[0][grid[0].length-2]
				}
			} else {
				new_grid.result[i][grid[0].length-1] = grid[i][grid[0].length-1];
			}
			if (neumann[3]){
				//if source
				if (grid[i][0] == 1){
					//apply dirichlet
					new_grid.result[i][0] = 1
				//else (=not source)
				} else {
					//apply neumann
					new_grid.result[i][0] = new_grid.result[i][1]
				}
			} else {
				new_grid.result[i][0] = grid[i][0];
			}
		}
		if (new_grid.error < max_error) break;
	}
	return new_grid;
}

function patch_grid(grid, neumann, flow=false, max_sim_error=0.1, max_loop=100){
	let new_grid = equilibrate(grid, neumann, flow, max_sim_error, max_loop);
	let frame=Float32Array.of(0,0,0,0);
	for (let i=0;i<grid.length;i++){
		frame[0]+=new_grid.result[i][grid[0].length-2];
		frame[2]+=new_grid.result[i][1];
	}
	for (let j=0;j<grid[0].length;j++){
		frame[1]+=new_grid.result[1][j];
		frame[3]+=new_grid.result[grid.length-2][j];
	}
	return {...new_grid, 'frame':frame};
}

function full_window(heat, flow, max_error=0.0001, max_loop=100){
	var new_flow=equilibrate(flow, [1,1,1,1], false, max_error, max_loop);
	//boundary conditions
	flow=new_flow.result;
	for (let i=0;i<flow.length;i++){
		flow[i][0]               =flow[i][1];
		//all edges are initially 0-valued
		flow[i][flow[0].length-1]=1;
	}
	for (let j=0;j<flow[0].length;j++){
		flow[0][j]            =flow[1][j];
		flow[flow.length-1][j]=flow[flow.length-2][j];
	}
	var new_heat= equilibrate(heat, [1,1,1,1], flow, max_error, max_loop);
	heat=new_heat.result;
	for (let i=0;i<heat.length;i++){
		heat[i][0]               =heat[i][0]==1? 1:heat[i][1];
		heat[i][heat[0].length-1]=heat[i][heat[0].length-2];
	}/*
	for (let j=0;j<heat[0].length;j++){
		heat[0][j]            =heat[1][j];
		heat[heat.length-1][j]=heat[heat.length-2][j];
	}*/
	return {'heat':heat, 'flow':flow}
}

function vec_add(a,b){
	return [a[0]+b[0],a[1]+b[1]];
}

function follow_border(flow, heat){
	var border={};
	var best=[];
	if (flow[1][1] == 0){
		border = {"met":[1, flow[1].lastIndexOf(0)], 'norm':[0,1]};
		best = [-10000,[1,border.met[1]+1]];
	} else {
		for(let i=1;i<flow.length-1;i++){
			if (flow[i][1]==0){
				border = {"met":[i,1], 'norm':[-1,0]};
				best=[-10000, [i-1,1]]
				break;
			}
		}
	}
	let candidates=[];
	total_grad = 0;
	for(let i=0;i<flow.length*flow[0].length;i++){

		//get values
		let air=vec_add(border.met,border.norm);
		let met=border.met;
		//ensure sanity
		if(air[0] > heat.length-2
						||
			met[0] > heat.length-2
						||
			air[0] <= 0
						||
			air[1] <= 0
						||
			met[0] <= 0 
						||
			met[1] <= 0
			) break;

		//update candidates /w sane values
		if (!candidates.length || (
			air[0] != candidates[candidates.length-1][1][0]
			       ||
			air[1] != candidates[candidates.length-1][1][1]
			)){
			//circular neighborhood
			var neighbors = [
				{'i':air[0],'j':air[1]-1},
				{'i':air[0]-1,'j':air[1]},
				{'i':air[0]+1,'j':air[1]},
				{'i':air[0],'j':air[1]+1}
			];
			let bricks = 0;
			let count = 0;
			let exhaust = 0;
			let exposed = 0
			for (let a=0; a<4; a++){
				let o=neighbors[a];
				if (flow[o.i][o.j] == 0){
					bricks += heat[o.i][o.j];
					count += 1;
				} else {
					exhaust += heat[o.i][o.j];
					exposed += 0.5;
				}
			}
			//TRUTH!
			let grad = bricks - heat[air[0]][air[1]];
			total_grad += grad;
			//GUESS!
			let estimate = exposed * bricks/count - exhaust;
			candidates.push([estimate - grad, air]);
			//guess is bad if a metal-edge would be encountered!
			let bad_guess = false;
			if ((estimate - grad)>best[0]){
				for (let i=-1;i<2;i+=2){
					for (let j=-1;j<2;j+=2){
						if(flow[air[0]+i][air[1]+j] == 0
								&&
							flow[air[0]+i][air[1]] //!= 0 resolves as true
								&&
							flow[air[0]][air[1]+j])
						{
							bad_guess=true;
							candidates[candidates.length-1].push(bad_guess);
							break;
						}
					}
				}
				if (! bad_guess){
					best = [estimate - grad, air];
				}
			}
		}

		//step values for next loop
		//rotate surface-normal vector
		let sideways=[border.norm[1],-border.norm[0]];
		//guess1
		g1=vec_add(met, sideways);
		//guess2
		g2=vec_add(air, sideways);

		//retrieve
		ismet1 = flow[g1[0]][g1[1]]==0;
		ismet2 = flow[g2[0]][g2[1]]==0;
		if (!ismet1){
			//!stay!
			//border.met = border.met
			border.norm = sideways;
			air=g1;
		} else if (ismet1 && ismet2){
			met=g2;
			border.met = g2;
			border.norm = [-sideways[0],-sideways[1]];
		} else if (ismet1 && !ismet2){
			met=g1;
			border.met = g1;
			air=g2;
			//border.norm = border.norm
		} else {//???
			throw(up);
		}

		//rerun loop with nice values
	}
	return {'total':total_grad, 'max':best[0], 'candidate':best[1]};
}

function extend_window(old_win, old, new_win, full){
	//helper-var: diff in start.x
	let start_diff=old_win.endpoints[0][0]-new_win.endpoints[0][0];
	var win = get_fresh_array(...new_win.dim);
	//*\\ off by 1 errors be lurking... //*\\
	//Every window comes with a FRAME we want to DISCARD!

	//for width of new window
	for (i=0;i<new_win.dim[0];i++){
		//if left of old window
		if(i<start_diff+1 ||
			//if right of old window (i is 0-indexed)
			i>(start_diff + old_win.dim[0] -2) ){
			//                              startpoint          + offset
			win[i]=Float32Array.of(...full[new_win.endpoints[0][0] + i].slice(
					//1st endpoint inclusive, 2nd exclusive. slice excludes last
					new_win.endpoints[0][1],new_win.endpoints[1][1] // this is fine.
				))
		//if INSIDE window
		} else {
			//selection as before
			win[i]=Float32Array.of(...full[new_win.endpoints[0][0]+i].slice(
					//excluding start of old_win is desirable now!
					new_win.endpoints[0][1],old_win.endpoints[0][1]+1
				//translate i into old-win x coordinate
				), ...old[i-start_diff].slice(1,old[0].length-1),
				//selection as before
				...full[new_win.endpoints[0][0]+i].slice(
					//exclude end of old_win, include end of new_win
					old_win.endpoints[1][1]-1,new_win.endpoints[1][1]
				)
			)
		}
	}
	return win;
}

function build_window(center, radius, limit){
	let x1 = max(center[0]-radius, 0);
	let y1 = max(center[1]-radius, 0);
	let x2 = min(center[0]+radius+1, limit[0]);
	let y2 = min(center[1]+radius+1, limit[1]);
	return {'dim':[x2-x1, y2-y1], 'endpoints':[[x1,y1],[x2,y2]]};
}

function extract(values, endpoints){
	return values.slice(endpoints[0][0], endpoints[1][0]).map((column) => column.slice(
			endpoints[0][1],endpoints[1][1]
		));
}

function dynamic_equilibrate(center, grid, flow=false, frame_error=0.001){
	let r=5;
	let new_win = build_window(center, r, [grid.length, grid[0].length]);
	let patch = extract(grid, new_win.endpoints)

	while(new_win.endpoints[0][0]!=0||new_win.endpoints[0][1]!=0||new_win.endpoints[1][0]!=grid.length||new_win.endpoints[1][1]!=grid[0].length){
		var old_win=new_win;
		//equilibrate window
		let neumann = [
			old_win.endpoints[1][0]>=grid.length-1,// && !flow,
			old_win.endpoints[1][1]>=grid[0].length-1,
			old_win.endpoints[0][0]<=0,
			old_win.endpoints[0][1]<=0 && !flow
		];
		flow_window = flow? extract(flow, old_win.endpoints) : false;
		var equilibrated = patch_grid(patch, neumann, flow_window, 0.000001, 1000);
		if(flow){
			update_table(equilibrated.result, old_win.endpoints);
		}
		//evaluate window
		//directional error specs?
		if (equilibrated.frame.every((diff) => diff < frame_error)){
			break;
		}

		//primitive. Could extend window in specific direction, based on frame-error?
		r += ~~(r/2);
		new_win = build_window(center, r, [grid.length, grid[0].length]);
		//patch = extract(grid, new_win.endpoints)
		patch = extend_window(old_win, equilibrated.result, new_win, grid);
	}
	return {
			'endpoints':old_win.endpoints,
			...equilibrated
		};
}

function get_base_T(flow, heat){
	let C = (c_metal+c_air+c_contact*cell/2)/(cell*L * c_metal*c_air*2*L*c_contact*cell*L);
	let border = get_surface_T(flow, heat);
	let T = Q_in/border.total/C;
	//turn to metal!
	flow[border.candidate[0]][border.candidate[1]] = 0;
}

function report_T(total_grad){
	let C = (c_metal+c_air+c_contact*cell/2)/(cell*L * c_metal*c_air*2*L*c_contact*cell*L);
	//Q_in is user reported!
	let T = Q_in/C / total_grad + room_T;
	document.getElementByID("reportT").text=""+T+" °C";
}

function get_table_elm(grid){
    var table = document.createElement("TABLE");
        
    for(let j = 0; j < grid[0].length; j++) {
        var row = table.insertRow(j);
        for(let i=0; i<grid.length; i++){
        	let val = grid[i][grid[0].length-j-1]
        	let cell =row.insertCell(i)
        	cell.bgColor = val==0? "#000" : "#"+ (0x1000000 + 0x10000*~~(0xff* val) + ~~(0xff * (1-val))).toString(16).slice(1,7);
        	cell.title=val;
        }
    }
    return table;
}

function update_table(grid, endpoints){
	var table = document.getElementById('result').children[0].children[0].children
	for(let j=0; j<height+1;j++){
		if (j<endpoints[0][1] || j>endpoints[1][1]){
			table[height-j].classList.remove('marked');
		} else {
			table[height-j].classList.add('marked');
			for (let i=0; i<length+2;i++){
				if (i<endpoints[0][0] || i>endpoints[1][0]){
					table[height-j].children[i].classList.remove('marked');
				} else {
					table[height-j].children[i].classList.add('marked');
					val = grid[i-endpoints[0][0]][j-endpoints[0][1]]
					//table[height-j].children[i].bgColor = "#"+ (0x1000000 + 0x10000*~~(0xff* val) + ~~(0xff * (1-val))).toString(16).slice(1,7);
					table[height-j].children[i].bgColor = val==0? "#000" : "#"+ (0x1000000 + 0x10000*~~(0xff* val) + ~~(0xff * (1-val))).toString(16).slice(1,7);
        			table[height-j].children[i].title=val;
				}
			}
		}
	}
}

function RunSimulation(){
	let repeat = ~~document.getElementById('nofsteps').value;
	for (let i=0;i<repeat;i++){
		setTimeout(function(){
			res = full_window(GlobalState.heat, GlobalState.flow, 0.000001,10000);
			let heat=res.heat;
			let flow=res.flow;
			let border=follow_border(flow, heat);
			candidate=border.candidate;
			flow[candidate[0]][candidate[1]] = 0;
			GlobalState = full_window(heat, flow, 0.000001,10000);
			Draw();
		},i*200);
	}
}

function run_optimisation(width, height){
	function _backpatch(values, patch){
		let left= patch.endpoints[0][0]>0 ? 1:0;
		let right=patch.result.length - (patch.endpoints[1][0]<width? 1:0);
		let bottom=patch.endpoints[0][1]>0 ? 1:0
		let top = patch.endpoints[1][1]<height? 1:0
		for (let i=left; i<right; i++){
			values[patch.endpoints[0][0]+i] = [
				...values[left+i].slice(0,patch.endpoints[0][1]+bottom),
				...patch.result[i].slice(bottom,patch.result[0].length-top),
				...values[left+i].slice(patch.endpoints[1][1]-top, values[0].length)
				]
		}
		return values;
	}
	var flow = prepare_air(width, height, [Math.floor(4*width/10),Math.floor(6*width/10)]);
	var heat = prepare_fire(width, height, [Math.floor(4*width/10),Math.floor(6*width/10)]);
	let win = full_window(heat, flow, 0.000000001, 10000);
	document.getElementById("result").innerHTML = get_table_elm(win.heat).outerHTML;
	heat=win.heat;
	flow=win.flow;
	let last_total = 0;
	//already metal
	let candidate = [15,1];

	for(let count=0;count<600;count++){

		//turn to metal!
		flow[candidate[0]][candidate[1]] = 0;
		console.log(candidate)
/*
		let flow_patch = dynamic_equilibrate(candidate, flow);
		flow = _backpatch(flow, flow_patch);
		let heat_patch = dynamic_equilibrate(candidate, heat, flow=flow);
		heat = _backpatch(heat, heat_patch);
*/
		let res = full_window(heat, flow, 0.00001,10000);
		heat=res.heat
		flow=res.flow
		let border = follow_border(flow,heat);
		last_total = border.total;
		candidate=border.candidate;
	}
	return{'heat':heat, 'flow':flow}
	return full_window(heat, flow);
}

//document.getElementById("result").innerHTML = get_table_elm(run_optimisation(width, height).heat).outerHTML;
/*document.getElementById("result").innerHTML = get_table_elm(
	full_window(
		prepare_fire(width, height, [Math.floor(width/4), Math.floor(3*width/4)]),
		prepare_air(width, height, [Math.floor(width/4), Math.floor(3*width/4)]),
		0.001,
		1000
	).heat).outerHTML;*/

function prepareSimulation(){
	cell = ~~document.getElementById("cell").value;
	var width= Math.floor(~~document.getElementById("width").value / cell);
	var height=Math.floor(~~document.getElementById("height").value / cell);
	var source1=Math.floor(~~document.getElementById("clearence").value / cell);
	var source2=Math.floor(~~document.getElementById("source_width").value / cell);
	var source=[source1, source2+source1];

	// recalculate everything
   metal_to_air = c_air * c_contact * cell*L /(c_metal + c_air + c_contact * cell/2)
   air_to_metal = c_metal * c_contact * cell*L /(c_metal + c_air + c_contact * cell/2)
	extra_air = v * 2 * cell**2 * 1.01 * 1.15 / (L*c_air);

	return {'flow':prepare_air(width, height, source), 'heat':prepare_fire(width, height, source)};
}

function ResetSimulation(){
	GlobalState = prepareSimulation();
	Draw();
}

ResetSimulation();

function Draw(){
	let selection=document.querySelector("input[type='radio']:checked").value
	document.getElementById("result").innerHTML = get_table_elm(GlobalState[selection]).outerHTML
}