function r = vres(test, quantity, value, ref, tol)
%VRES  One validation record: relative error and pass/fail against tolerance.
r.test = test; r.quantity = quantity; r.value = value; r.ref = ref;
r.err = abs(value - ref) / max(abs(ref), eps);
r.tol = tol; r.pass = r.err <= tol;
end
