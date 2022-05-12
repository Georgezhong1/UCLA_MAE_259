function [F,Jacob] = Grad_Hs_EE_rod(q_new,EA,EI,GJ,voronoiRefLen,refLen, ...
    kappa,a1_new,a2_new,tangent_new,delta_mk,m1,m2,g,m,nv)
    F = zeros(4*nv-1,1);
    Jacob = zeros(4*nv-1,4*nv-1);

    for i = 1:nv-1
        ind = [4*i-3:4*i-1,4*i+1:4*i+3]';
        q_new(ind(1:3));
        q_new(ind(4:6));

        [dF, dJ] = gradEs_hessEs(q_new(4*i-3:4*i-1)',q_new(4*i+1:4*i+3)',refLen,EA);

        F(ind)= F(ind)+dF;
        Jacob(ind,ind) = Jacob(ind,ind)+dJ;
    end

    for i = 2:nv-1
        ind = (4*i-7 : 4*i+3);
        node0 = q_new(ind(1:3))';
        node1 = q_new(ind(5:7))';
        node2 = q_new(ind(9:11))';

        [dF,dJ] = gradEb_hessEb(node0,node1,node2, m1(i-1,:), m2(i-1,:), ...
            m1(i,:), m2(i,:), kappa(i-1,:), refLen, EI);
        F(ind) = F(ind)+dF;
        Jacob(ind,ind) = Jacob(ind,ind)+dJ;

        [dF,dJ] = gradEt_hessEt(node0,node1, ...
            node2,q_new(ind(4)),q_new(ind(8)),delta_mk(i-1),refLen,GJ);
        F(ind)= F(ind)+dF;
        Jacob(ind,ind) = Jacob(ind,ind)+dJ;
    end

    for i = 1:nv
        m_tmp = m(4*i-3:4*i-1,4*i-3:4*i-1);
        F(4*i-3:4*i-1) = F(4*i-3:4*i-1)- m_tmp*g';
    end

end
    