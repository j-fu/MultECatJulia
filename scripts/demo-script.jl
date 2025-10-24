module DemoScript
using MultECatJulia
using MySubPackage

function main()
    MultECatJulia.greet()
    return MySubPackage.greet()
end

end
