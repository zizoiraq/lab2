function parseCoefficients(str) {
    return str.split(',').map(function(x){return parseFloat(x.trim());});
}

function checkStability() {
    const den = parseCoefficients(document.getElementById('den').value);
    const roots = math.roots(den);
    const stable = roots.every(r => math.re(r) < 0);
    const out = stable ? 'System is stable' : 'System is unstable';
    document.getElementById('output').innerText = out + '\nRoots:\n' + JSON.stringify(roots, null, 2);
}

function computeRouth(coeffs) {
    let n = coeffs.length - 1;
    let rows = [];
    rows.push(coeffs.filter((_, i) => i%2===0));
    rows.push(coeffs.filter((_, i) => i%2===1));
    if (rows[1].length < rows[0].length) rows[1].push(0);
    for (let i=2; i<=n; i++) {
        rows[i] = [];
        for (let j=0; j<rows[i-1].length-1; j++) {
            let a = rows[i-1][0];
            if (a === 0) a = 1e-6; // avoid div by zero
            let b = rows[i-1][j+1];
            let c = rows[i-2][0];
            let d = rows[i-2][j+1] || 0;
            rows[i][j] = ((a*d) - (c*b)) / a;
        }
        if (rows[i].length < rows[i-1].length) rows[i].push(0);
    }
    return rows;
}

function showRouth() {
    const den = parseCoefficients(document.getElementById('den').value);
    const routh = computeRouth(den);
    let text = 'Routh Table:\n';
    routh.forEach((row, i) => {
        text += 'Row ' + i + ': ' + row.map(x=>x.toFixed(3)).join(' ') + '\n';
    });
    document.getElementById('output').innerText = text;
}

function plotRootLocus() {
    const num = parseCoefficients(document.getElementById('num').value);
    const den = parseCoefficients(document.getElementById('den').value);
    let re = [];
    let im = [];
    for (let k=0; k<=100; k+=0.5) {
        let poly = den.map((d,i)=>d).slice();
        for (let i=0;i<num.length;i++) {
            poly[poly.length-num.length+i] += k*num[i];
        }
        let roots = math.roots(poly);
        roots.forEach(r => {re.push(math.re(r)); im.push(math.im(r));});
    }
    const data = [{x: re, y: im, mode:'markers', type:'scatter', name:'roots'}];
    const layout = {xaxis:{title:'Real'}, yaxis:{title:'Imaginary'}, title:'Root Locus'};
    Plotly.newPlot('plot', data, layout);
    document.getElementById('output').innerText = '';
}

function plotBode() {
    const num = parseCoefficients(document.getElementById('num').value);
    const den = parseCoefficients(document.getElementById('den').value);
    let w = [];
    for (let i=-2;i<=2;i+=0.02) {
        w.push(Math.pow(10,i));
    }
    let mag = [];
    let phase = [];
    w.forEach(freq => {
        const s = math.complex(0,freq);
        const Ns = num.reduce((acc,c,i)=>math.add(acc, math.multiply(c, math.pow(s, num.length-1-i))), math.complex(0,0));
        const Ds = den.reduce((acc,c,i)=>math.add(acc, math.multiply(c, math.pow(s, den.length-1-i))), math.complex(0,0));
        const G = math.divide(Ns, Ds);
        mag.push(20*Math.log10(math.abs(G)));
        phase.push(math.atan2(math.im(G), math.re(G))*180/Math.PI);
    });
    const data = [
        {x:w, y:mag, type:'scatter', name:'Magnitude'},
        {x:w, y:phase, type:'scatter', yaxis:'y2', name:'Phase'}
    ];
    const layout = {
        title:'Bode Plot',
        xaxis:{type:'log',title:'Frequency (rad/s)'},
        yaxis:{title:'Magnitude (dB)'},
        yaxis2:{title:'Phase (deg)', overlaying:'y', side:'right'}
    };
    Plotly.newPlot('plot', data, layout);
    document.getElementById('output').innerText = '';
}

function pidResponse() {
    const num = parseCoefficients(document.getElementById('num').value);
    const den = parseCoefficients(document.getElementById('den').value);
    const [kp,ki,kd] = parseCoefficients(document.getElementById('pid').value);
    const order = den.length-1;
    // companion form
    let A = math.zeros(order, order);
    for (let i=0;i<order-1;i++) A.set([i, i+1], 1);
    for (let i=0;i<order;i++) A.set([order-1, i], -den[i]/den[order]);
    let B = math.zeros(order,1);
    B.set([order-1,0],1/den[order]);
    let C = math.matrix(num.length===order+1 ? num.slice(0,order) : Array(order).fill(0));
    C = math.divide(C, den[order]);
    let D = num.length===order+1 ? num[num.length-1]/den[order] : 0;
    const dt = 0.01;
    const steps = 1000;
    // discretize
    let Ad = math.add(math.identity(order), math.multiply(A, dt));
    let Bd = math.multiply(B, dt);
    let x = math.zeros(order,1);
    let y = 0;
    let tdata=[]; let ydata=[];
    let integral = 0; let prevErr = 0;
    for (let i=0;i<steps;i++) {
        let t = i*dt;
        let err = 1 - y;
        integral += err*dt;
        let derivative = (err - prevErr)/dt;
        prevErr = err;
        let u = kp*err + ki*integral + kd*derivative;
        x = math.add(math.multiply(Ad, x), math.multiply(Bd, u));
        y = math.add(math.multiply(C, x).toArray().reduce((a,b)=>a+b,0), D*u);
        tdata.push(t); ydata.push(y);
    }
    const data = [{x:tdata, y:ydata, type:'scatter', name:'y(t)'}];
    const layout={title:'PID Step Response', xaxis:{title:'Time (s)'}, yaxis:{title:'Output'}};
    Plotly.newPlot('plot', data, layout);
    document.getElementById('output').innerText = '';
}

function steadyStateError() {
    const num = parseCoefficients(document.getElementById('num').value);
    const den = parseCoefficients(document.getElementById('den').value);
    const N0 = num[num.length-1];
    const D0 = den[den.length-1];
    const Kp = N0/D0;
    const ess = 1/(1+Kp);
    document.getElementById('output').innerText = 'Steady state error for step input: ' + ess.toFixed(4);
}
