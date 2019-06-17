/********************************************************************
geometry.js
********************************************************************/
function angles_perpendicular(angle1, angle2){ //6
    /* Returns True if angles are perpendicular or almost perpendicular.
       There is a margin of +-0.1 radians in which they are still
       considered perpendicular. */
    if(Math.abs(angle2 - angle1 - (3.1416 / 2)) < 0.1 || Math.abs(angle2 - angle1 + (3.1416 / 2)) < 0.1){
    	return true;
    }else{
    	return false;
    }   
}

function intersection(hline, vline){ //7
    /*Returns the intersection point of a (nearly) horizontal line
       with a (nearly) vertical line. Results may be wrong in
       other cases.*/     
    rho1 = hline[0];
    theta1 = hline[1];
    rho2 = vline[0];
    theta2 = vline[1];
    y_coor = (rho1 * Math.cos(theta2) - rho2 * Math.cos(theta1)) / Math.sin(theta1 - theta2);
    x_coor = (rho2 - y_coor * Math.sin(theta2)) / Math.cos(theta2);
    point = {x:Math.round(x_coor), y:Math.round(y_coor)};
    return point;
}

function distance(p1,p2){
	return Math.sqrt((p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y));
}

function closer_points_rel(p1,p2,offset_ratio,abs_offset){
	let dx = p2.x - p1.x;
	let dy = p2.y - p1.y;
	let xoff = abs_offset*dx / (Math.abs(dx) + Math.abs(dy));
    let yoff = abs_offset*dy / (Math.abs(dx) + Math.abs(dy));
    let k = (1 - offset_ratio)/2;
    return new Array({x:Math.round(p1.x+(dx*k)+xoff), y:Math.round(p1.y+(dy*k)+yoff)}, {x:Math.round(p2.x-(dx*k)-xoff), y:Math.round(p2.y-(dy*k)-yoff)});
}

function mean_corn(a,b,c,d){
	return ((a+b+c+d)/4);
}
/********************************************************************/
