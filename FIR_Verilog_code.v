module FIR_Lowpass_Filter (
    input wire clk,
    input wire [15:0] data_in,
    output reg [31:0] data_out
);

    // ===== Type Definitions =====
    reg signed [15:0] delay_line [0:31];
    wire signed [26:0] products [0:31];
    reg signed [31:0] sum;
    
    parameter  N = 32; // Number o
    
    // ===== Coefficients =====
    wire signed [10:0] coeffs [0:31];
    assign coeffs[0] = -17;
    assign coeffs[1] = -20;
    assign coeffs[2] = -26;
    assign coeffs[3] = -31;
    assign coeffs[4] = -29;
    assign coeffs[5] = -15;
    assign coeffs[6] = 20;
    assign coeffs[7] = 82;
    assign coeffs[8] = 174;
    assign coeffs[9] = 294;
    assign coeffs[10] = 437;
    assign coeffs[11] = 591;
    assign coeffs[12] = 741;
    assign coeffs[13] = 873;
    assign coeffs[14] = 971;
    assign coeffs[15] = 1023;
    assign coeffs[16] = 1023;
    assign coeffs[17] = 971;
    assign coeffs[18] = 873;
    assign coeffs[19] = 741;
    assign coeffs[20] = 591;
    assign coeffs[21] = 437;
    assign coeffs[22] = 294;
    assign coeffs[23] = 174;
    assign coeffs[24] = 82;
    assign coeffs[25] = 20;
    assign coeffs[26] = -15;
    assign coeffs[27] = -29;
    assign coeffs[28] = -31;
    assign coeffs[29] = -26;
    assign coeffs[30] = -20;
    assign coeffs[31] = -17;

    // ===== Multiply Stage =====
    genvar i;
    
    generate
        for (i = 0; i < 32 ; i = i + 1) begin : mult_stage
            assign products[i] = delay_line[i] * coeffs[i];
        end
    endgenerate

    // ===== Main Logic =====
    integer j;
    always @(posedge clk) begin
        // Shift delay line
        for (j = 31; j > 0; j = j - 1) begin
            delay_line[j] <= delay_line[j - 1];
        end
        delay_line[0] <= data_in;

        // Single-stage sum
        sum = 0;
        for (j = 0; j < 32; j = j + 1) begin
            sum = sum + products[j];
        end

        // Assign to output
        data_out <= sum;
    end

endmodule
