from pathlib import Path
import textwrap
import numpy as np
import pytest
import mkin4py


@pytest.mark.parametrize('method', ['dense', 'qmr'])
def test_original_readme_example(method):
    readme = (Path(__file__).parents[1] / 'README.md').read_text()
    code = textwrap.dedent(readme[readme.index('        import mkin4py'):readme.index('  - **Evaluation**:')])
    exec(compile(code, 'README.md', 'exec'), {})
    mkin4py.mkmodel.reset_model(seed=0)
    result = mkin4py.solver.solve.rk4(linear_solver=method, max_restarts=3)
    assert result['success'], result
    expected = [.439342950,.00119743307,.107992516,.110447591,2.99730332e-9,
                7.70711567e-10,.100256049,.0978269843,.132727193,.00964368428,
                3.28419897e-8,4.59118359e-13,.000565562897,1.12150752e-15]
    np.testing.assert_allclose(result['coverage'][6:,0], expected, rtol=3e-6, atol=1e-12)
    np.testing.assert_allclose(result['rates'][:6,0],
        [-42.425219,-24.9532633,12.9732696,.000558723011,23.9588699,23.9588699], rtol=3e-6)
    np.testing.assert_allclose(result['rates'][6:], 0, atol=1e-8)
