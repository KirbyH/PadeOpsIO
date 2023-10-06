import os 
import sys
import numpy as np

# sys.path.append(r'./padeops_io/')  # why is this our cwd? 
import padeopsio as pio

# testIO = BudgetIO('some_dir', verbose=False)
testIO = pio.BudgetIO(r'./test_data', verbose=True)

pio.PlotsIO.xy_slice(testIO, ['ubar', 'uwake'], z=5)

budgets = testIO.existing_budgets()
testIO.existing_terms(budgets[0])
testIO.read_budgets(budget_terms=['ubar'])
testIO.read_budgets(budget_terms=['ubar', 'vbar', 'wbar', (0, 1), '(0, 1)', 'pbar'])
testIO.calc_wake()
# testIO.write_npz(budget_terms='all', filename='overwrite_test', overwrite=True)

InflowParser.inflow_budgets(testIO)


print(testIO.__dict__)
