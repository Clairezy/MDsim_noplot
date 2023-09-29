//
//  MDsim.swift
//  MDsim_noplot
//
//  Created by Claire on 9/19/23.
//
import Foundation

// initialize global variables
var L = 0                 // chain size
let Natom = 36            // # of atoms
let Nmax = 200            // max # of atoms (for arrays)
var x = [Double](repeating: 0.0, count: Nmax) // array to store positions
var fx = [[Double]](repeating: [Double](repeating: 0.0, count: 2), count: Nmax) // array to store forces

// maaxwellian velocity distribution
func maxwellianVelocity(T: Double) -> Double {
    var s = 0.0
    // sum 24 random numbers between 0 and 1
    for _ in 0..<36 {
        s += Double.random(in: 0.0..<1.0)
    }
    // adjust and scale to get Maxwellian fixed temp velocity
    return (s - 18.0) * sqrt(T)
}

// main simulation function
func main() {
    var t1 = 0, t2 = 1, t = 0 // time variables
    let h = 0.00002           // most stable time step size
    let hover2 = h / 2.0      // half of the time step
    let Nstep = 5000          // # of simulation steps
    let Nprint = 100          // print results every Nprint steps
    let Ndim = 1              // # of dimensions
    let Tinit = 150.0         // initial temperature
    var vx = [Double](repeating: 0.0, count: Nmax) // array to store velocities

    // calculate the simulation box size
    L = Int(pow(Double(Natom), 1.0 / Double(Ndim)))

    print("Natom = \(Natom) L = \(L)")
    
    // initialize particle positions and velocities
    for ix in 0..<L {
        x[ix] = Double(ix)
        vx[ix] = maxwellianVelocity(T: Tinit)
    }

    var KE = 0.0 // kinetic energy
    var PE = 0.0 // potential energy

    // calculate initial potential energy using Forces function
    PE = Forces(t1, PE)
    for i in 0..<Natom {
        KE += (vx[i] * vx[i]) / 2.0 // calculate initial kinetic energy
    }

    // calculate initial scaled temperature
    let Tinit_scaled = 2.0 * KE / (3.0 * Double(Natom))
    print("Initial Temperature Scaled (Tinit) = \(Tinit_scaled)")
    
    print("\(t) PE = \(PE) KE = \(KE) PE + KE = \(PE + KE)")

    // arrays to store kinetic and potential energies over time
    var kineticEnergies = [Double]()
    var potentialEnergies = [Double]()

    // main simulation loop
    for t in 1..<Nstep {
        // update particle positions and calculate potential energy
        for i in 0..<Natom {
            PE = Forces(t1, PE)
            x[i] = x[i] + h * (vx[i] + hover2 * fx[i][t1])
            if x[i] <= 0.0 {
                x[i] = x[i] + Double(L)
            }
            if x[i] >= Double(L) {
                x[i] = x[i] - Double(L)
            }
        }

        // update potential energy and kinetic energy
        PE = Forces(t2, PE)
        KE = 0.0
        for i in 0..<Natom {
            vx[i] = vx[i] + hover2 * (fx[i][t1] + fx[i][t2])
            KE += (vx[i] * vx[i]) / 2.0
        }

        // print energies at regular intervals
        if t % Nprint == 0 {
            print("\(t) PE = \(PE) KE = \(KE) PE + KE = \(PE + KE)")
        }

        // store energies for later
        kineticEnergies.append(KE)
        potentialEnergies.append(PE)

        // swap time indices t1 and t2
        (t1, t2) = (t2, t1)
    }
    
    // calculate final temperature
    let T = 2.0 * KE / (3.0 * Double(Natom))
    print("Final Temperature (T) = \(T)")
    
    // evaluate time-averaged energies at equilibrium
    let numStepsForAveraging = 1000 // # of last steps for averaging
    var timeAveragedKE = 0.0
    var timeAveragedPE = 0.0

    // loop to evaluate time-averaged energies
    for t in (Nstep - numStepsForAveraging)..<(Nstep-1) {
            timeAveragedKE += kineticEnergies[t]
            timeAveragedPE += potentialEnergies[t]
    }

    // calculate and print time-averaged energies
    timeAveragedKE /= Double(numStepsForAveraging)
    timeAveragedPE /= Double(numStepsForAveraging)
    print("Time-Averaged Kinetic Energy = \(timeAveragedKE)")
    print("Time-Averaged Potential Energy = \(timeAveragedPE)")
}

// function to calculate forces between particles and update potential energy
func Forces(_ t: Int, _ PE: Double) -> Double {
    let r2cut = 9.0 // cutoff radius squared for force calculations
    var PE = 0.0 // initialize potential energy

    // initialize force arrays to zero
    for i in 0..<Natom {
        fx[i][t] = 0.0
    }

    // calculate forces and potential energy between particles
    for i in 0..<(Natom - 1) {
        for j in (i + 1)..<Natom {
            var dx = x[i] - x[j]

            // apply periodic boundary conditions
            if abs(dx) > 0.5 * Double(L) {
                dx = dx < 0 ? dx + Double(L) : dx - Double(L)
            }

            let r2 = dx * dx // distance squared between particles

            // check if particles are within the cutoff radius
            if r2 < r2cut {
                if r2 == 0.0 {
                    dx = 0.0001
                }
                let invr2 = 1.0 / r2
                let fijx = 48.0 * (pow(invr2, 3) - 0.5) * pow(invr2, 3) * invr2 * dx
                fx[i][t] = fx[i][t] + fijx
                fx[j][t] = fx[j][t] - fijx
                PE = PE + 4.0 * pow(invr2, 3) * (pow(invr2, 3) - 1.0)
            }
        }
    }

    return PE
}
