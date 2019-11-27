function getelectricfieldproperty(e::ElectricField{1},sym::Symbol)
    n,m = size(getfield(e,:values))
    if sym==:left
        start_ind = getfield(e,:start_inds)[1]
        values = getfield(e,:values)
        inds = [0n÷3+start_ind,1n÷3+start_ind,2n÷3+start_ind]
        p = values[inds,:]
    elseif sym==:right
        stop_ind = getfield(e,:stop_inds)[1]
        values = getfield(e,:values)
        inds = [0n÷3+stop_ind,1n÷3+stop_ind,2n÷3+stop_ind]
        p = values[inds,:]
    else
        return nothing
    end
end
