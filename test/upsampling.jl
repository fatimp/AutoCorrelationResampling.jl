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

function contract(array, times)
    slices = (array[:,:,k] for k in 1:times:size(array, 3))
    return reduce((x, y) -> cat(x, y; dims = (3)), slices)
end

@testset "Test upsampling (x3 .. x6)" begin
    for n in 3:6
        cut_data = contract(data, n)
        ac = autocorrelation(cut_data)
        filter = filter_coeffs(n)
        ac_upsampled = ac_upsample(ac, filter)
        ac_upsampled_ft = fft(ac_upsampled)

        @test (ac_upsampled |> imag |> maximum) < +1e-8
        @test (ac_upsampled |> real |> minimum) > -1e-8
    end
end
