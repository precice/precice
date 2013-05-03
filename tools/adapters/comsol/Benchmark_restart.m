% make this state as initial state
init = asseminit(fem,'init',fem0.sol,'solnum',solsize(fem0.sol));
fem.sol = init; %u = asseminit(fem,'init',fem0.sol,'solnum',solsize(fem0.sol));