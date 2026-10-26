"""Independent small-array checks for the feasibility mathematics."""
import json
import unittest
import numpy as np
import study as s

class GramTests(unittest.TestCase):
    def test_projector_gram(self):
        rng=np.random.default_rng(1205);u=rng.normal(size=(3,4));v=rng.normal(size=(2,4))
        p=np.einsum('bi,bj->bij',u,u).reshape(3,-1)
        q=np.einsum('bi,bj->bij',v,v).reshape(2,-1)
        np.testing.assert_allclose(s.gram(u,v),p@q.T,atol=1e-12)

    def test_variable_projection_gradient(self):
        rng=np.random.default_rng(3107)
        u0=rng.normal(size=(1,3,4));u=rng.normal(size=(1,2,4));w0=rng.normal(size=(1,3,2))
        w=s.weights_for(u,u0,w0);g=s.gradient(u,w,u0,w0)
        step=1e-6;fd=np.empty_like(u)
        for index in np.ndindex(u.shape):
            plus=u.copy();minus=u.copy();plus[index]+=step;minus[index]-=step
            fp=s.factor_loss(plus,s.weights_for(plus,u0,w0),u0,w0).sum()
            fm=s.factor_loss(minus,s.weights_for(minus,u0,w0),u0,w0).sum()
            fd[index]=(fp-fm)/(2*step)
        np.testing.assert_allclose(g,fd,atol=2e-6,rtol=2e-6)

    def test_coupled_removal_and_gradient(self):
        rng=np.random.default_rng(9173);x=rng.normal(size=(2,3,5));raw=rng.normal(size=6)
        e=rng.normal(size=(5,5));e=(e+e.T)/2
        cross=(x.reshape(6,5)@x.reshape(6,5).T).reshape(2,3,2,3)
        a=np.einsum('rcp,pq,rdq->rcd',x,e,x)
        value,g=s.coupled_value_gradient(raw,cross,a,1.)
        z=raw.reshape(2,3);z=z/np.linalg.norm(z,axis=1,keepdims=True)
        y=np.einsum('rc,rcp->rp',z,x)
        direct=np.sum((e-y.T@y)**2)-np.sum(e**2)
        self.assertAlmostEqual(value,direct,places=10)
        step=1e-6;fd=np.empty(6)
        for j in range(6):
            plus=raw.copy();minus=raw.copy();plus[j]+=step;minus[j]-=step
            fd[j]=(s.coupled_value_gradient(plus,cross,a,1.)[0]-
                   s.coupled_value_gradient(minus,cross,a,1.)[0])/(2*step)
        np.testing.assert_allclose(g,fd,atol=1e-6,rtol=2e-7)
        with self.assertRaises(ArithmeticError):s.coupled_value_gradient(np.zeros(6),cross,a,1.)

    def test_wl_gradient_fixture(self):
        case=json.loads((s.HERE/'symbolic.json').read_text())['gradient_fixture']
        args=[np.array(case[k]) for k in ['u','w','u0','w0']]
        np.testing.assert_allclose(s.gradient(*args),case['gradient'],atol=1e-12)

    def test_greedy_projection_and_edges(self):
        rng=np.random.default_rng(491);u0=rng.normal(size=(1,4,5));w0=rng.normal(size=(1,4,2))
        errors=[]
        for count in [1,2,3,4]:
            u,indices=s.greedy_dictionary(u0,w0,count)
            self.assertEqual(len(set(indices[0])),count)
            errors.append(float(s.factor_loss(u,s.weights_for(u,u0,w0),u0,w0).sum()))
        self.assertTrue(np.all(np.diff(errors)<=1e-9));self.assertAlmostEqual(errors[-1],0,places=9)
        for bad in [0,5,1.5,True]:
            with self.assertRaises(ValueError):s.greedy_dictionary(u0,w0,bad)

if __name__=='__main__':unittest.main()
