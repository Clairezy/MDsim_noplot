//
//  ContentView.swift
//  MDsim_noplot
//
//  Created by Claire on 9/19/23.
//

import SwiftUI

struct ContentView: View {
    var body: some View {
        
        Button("Run Simulation") {
            main() // Call main function
            
        }
    }
}


struct ContentView_Previews: PreviewProvider {
    static var previews: some View {
        ContentView()
    }
}
