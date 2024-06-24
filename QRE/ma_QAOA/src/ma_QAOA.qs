// This is a program for simulating the ma-QAOA circuit. 
namespace ma_qaoa{
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Measurement;

    // Cost Hamiltonian
    operation cost_unitary(qubits: Qubit[], gamma: Double[], edges: Int[], weights: Double[], n_qubits: Int): Unit{
        mutable gamma_count : Int = 0;
        
        // RZZ Gates
        for flat_index in edges{
            let (i,j) = DivRemI(flat_index,n_qubits);
            Rzz(2.0*gamma[gamma_count]*weights[gamma_count],qubits[i],qubits[j]);
            set gamma_count += 1;
        }
    }

    // Mixer Hamiltonian
    operation mixer_unitary(qubits: Qubit[], beta: Double[]) : Unit{
        
        mutable beta_count : Int = 0;

        for qubit in qubits{
            Rx(2.0 * beta[beta_count],qubit);
            set beta_count += 1
        }
    }

    // Function to create the QAOA circuit and return the measurement result as an integer.
    operation ma_circuit(NQubits: Int, layers: Int, gammas: Double[], betas: Double[], edges: Int[], weights: Double[]) : Int {

        use q = Qubit[NQubits]; 
        mutable integer_result = 0;
        let n_edges = Length(edges);
        
        // State Preparation - |+>
        ApplyToEachA(H,q);
        
        // Applying the layers  
        for layer in 0..layers-1{

            // Message($"Layers {layer}");

            cost_unitary(q, gammas[layer*n_edges..(layer+1)*n_edges-1], edges, weights, NQubits);
            mixer_unitary(q, betas[layer*NQubits..(layer+1)*NQubits-1]);
        }
        // Return the bitstring as an integer.
        return MeasureInteger(q);
    }
}
