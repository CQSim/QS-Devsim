
def region_solve(device, region, absolute_error, relative_error, maximum_iterations, equation=None):
  '''
    run solve until desired region meets relative_error and absolute_error criteria
    if equation == "", then check only the region convergence
  '''
  # run up to the maximum number of iterations
  for i in range(maximum_iterations):
    try:
      # this solve statement should almost always succeed
      info = solve(type="dc", absolute_error=1e30, relative_error=1e10, maximum_iterations=1, info=True)
    except:
      raise
    # get the last iteration
    last_iteration = info['iterations'][-1]
    # A StopIteration exception is thrown if the name is not found
    test_device = next(x for x in last_iteration['devices'] if x['name'] == device)
    test_region = next(x for x in test_device['regions'] if x['name'] == region)
    if equation:
      test_equation = next(x for x in test_region['equations'] if x['name'] == equation)
      converged = (test_equation['absolute_error'] < absolute_error) and (test_equation['relative_error'] < relative_error)
    else:
      converged = (test_region['absolute_error'] < absolute_error) and (test_region['relative_error'] < relative_error)
    if converged:
      break
  # if you prefer to raise an exception for nonconvergence, do so here
  if not converged:
    raise RuntimeError("Convergence failure")
  return info

