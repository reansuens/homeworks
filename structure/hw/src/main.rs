use std::f64::consts::PI;

fn main() {

    const BR: f64 = 150.0;          // mm^2
    const MX: f64 = 100.0;          // kN·m
    const SY: f64 = 50.0e3;         // N

    const Y: [f64; 10] = [
        750.0, 603.6, 250.0, -250.0, -603.6,
        -750.0, -603.6, -250.0, 250.0, 603.6
    ];

    let ixx = 2.0 * BR * (
        750.0_f64.powi(2) +
        2.0 * 603.6_f64.powi(2) +
        2.0 * 250.0_f64.powi(2)
    );

    let mut sigma_z: Vec<f64> = Vec::with_capacity(Y.len());

    for (_, y_val) in Y.iter().enumerate() {
        let sigma = (MX / ixx) * y_val * 1e6;
        sigma_z.push(sigma);
    }

    let c = -(SY * BR) / ixx;

    let mut qb: Vec<f64> = Vec::with_capacity(10);

    let qb910 = c * Y[8];
    let qb101 = qb910 + c * Y[9];
    let qb12  = qb101 + c * Y[0];
    qb.push(qb12);

    for (i, _) in Y.iter().enumerate().take(3).skip(1) {
        let prev = qb.last().unwrap();
        qb.push(prev + c * Y[i]);
    }

    // Symmetry propagation (remaining values)
    qb.push(qb[1]);
    qb.push(qb[0]);
    qb.push(qb101);
    qb.push(qb910);
    qb.push(0.0);
    qb.push(qb910);
    qb.push(qb101);

    
    let a34  = 0.5 * 500.0 * 250.0;
    let a23  = 62500.0 + (45.0/360.0) * PI * 500.0_f64.powi(2)
               - 0.5 * 250.0 * 353.6;
    let a910 = a23;
    let a12  = 0.5 * 250.0 * 353.6
               + (45.0/360.0) * PI * 500.0_f64.powi(2);
    let a101 = a12;
    let a    = 500.0 * 1000.0 + PI * 500.0_f64.powi(2);

    let qso = (
        SY * 250.0
        - 2.0 * (
            2.0 * qb910 * a910 +
            2.0 * qb101 * a101 +
            2.0 * qb12  * a12  +
            2.0 * qb[1] * a23  +
            2.0 * qb[2] * a34
        )
    ) / (-2.0 * a);

    let mut q: Vec<f64> = Vec::with_capacity(qb.len());

    for (_, val) in qb.iter().enumerate() {
        q.push((val - qso).abs());
    }

    let names = [
        "q12","q23","q34","q45","q56",
        "q67","q78","q89","q910","q101"
    ];

    println!("\n--- Shear Flow Results (Absolute Values) ---\n");

    for (i, value) in q.iter().enumerate() {
        println!("{:<6} = {:>12.6} N/mm^2", names[i], value);
    }
}

