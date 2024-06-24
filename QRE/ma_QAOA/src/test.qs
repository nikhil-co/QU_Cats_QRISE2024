
namespace testing{

    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Measurement;
    
    operation SayHelloQ(dic: Int[], n: Int) : Int {

        use q = Qubit[10]; 
        mutable integer_result = 0;
        let b = [1,2,3,45,5];
        // State Preparation - |+>
        ApplyToEachA(H,q);
        // Message($"The integer is: {dic[0]}");

        for index in dic {
            let (i,j) = DivRemI(index,n);
            Message($"The index is: {i} {j} and the flat index is: {index}");
            Message($"{b[1..3]}")
        }

        // Return the bitstring as an integer.
        return MeasureInteger(q);
    }
}