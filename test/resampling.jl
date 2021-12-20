function read_data(name, dims, side)
    s = Tuple(side for _ in 1:dims)
    array = zeros(UInt8, s)

    open(name) do io
        read!(io, array)
    end

    return BitArray(array)
end

const data = read_data("ceramic200", 3, 200)

@testset "Test filters (x3 .. x6)" begin
    for n in 3:6
        filter = filter_coeffs(n)
        @test 0.98 < AutoCorrelationResampling.csum(filter) < 1.02
    end
end

autocorrelation(array) = (array |> fft .|> abs2 |> ifft |> real) / length(array)

@testset "Test resampling" begin
    for n in [2, 3, 2//3, 3//2]
        ac = autocorrelation(data)
        ac_resampled = ac_resample(ac, n)
        ac_resampled_ft = fft(ac_resampled)

        @test (ac_resampled |> imag |> maximum) < +1e-8
        @test (ac_resampled |> real |> minimum) > -1e-8
    end
end
