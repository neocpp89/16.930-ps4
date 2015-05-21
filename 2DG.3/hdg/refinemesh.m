function [mout, mout1] = refinemesh(m, nref)
    mout = m;
    [mout.p, mout.t] = uniref(m.p, m.t, nref);
    fd=@(p) drectangle(p, 0, 1, 0, 1);
    [mout.f,mout.t2f] = mkt2f(mout.t);

    mout.f(mout.f(:, 4) <= 0, 4) = -1;

    mout.fcurved = (mout.f(:,4)<0);
    ic = find(mout.fcurved);
    mout.tcurved = repmat(false, size(mout.t,1),1);
    mout.tcurved(mout.f(ic,3)) = true;

    mout.porder = m.porder;
    [mout.plocal,mout.tlocal] = uniformlocalpnts(mout.porder);
    mout.dgnodes = createnodes(mout,fd);
end
