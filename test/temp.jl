module Reqs
    module META
        module METAMETA
            using Requires
            function __init__()
                @require SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b" @eval begin
                    hello(x,y,z) = "HI!"
                end
            end
            hello(x,y) = "HI!"
        export hello
    end
    using .METAMETA
    METAMETA.hello(x) = "HI!"
    export hello
end
using .META
META.hello() = "HI!"
export hello
end
