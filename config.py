import platform


RESULTS_PATH = None
print('Platform={}'.format(platform.node()))
if platform.node() == 'your_computer':
    REPO_PATH = '/code/SpectralDynamics'
    RESULTS_PATH = '/data'
    DB_PATH = '/data/=analysis/SpectralDynamics_target_v0.1.db'

elif platform.node() == 'different_computer':
    REPO_PATH = '/diff_code/SpectralDynamics'
    RESULTS_PATH = '/diff_data'
    DB_PATH = '/diff_data/=analysis/SpectralDynamics_target_v0.1.db'
