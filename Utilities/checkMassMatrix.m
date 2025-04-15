function checkMassMatrix(M)
    try chol(M);
    catch
            disp('Mass Matrix is not symmetric positive definite')
    end
end