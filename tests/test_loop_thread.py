import unittest
from unittest.mock import patch, MagicMock
from scripts.optimizer import LoopThread

class TestLoopThread(unittest.TestCase):
    def setUp(self):
        self.loop_thread = LoopThread()

    @patch("subprocess.Popen")
    def test_run_datcom(self, mock_popen):
        # Mock subprocess.Popen
        mock_process = MagicMock()
        mock_process.communicate.return_value = ("stdout", "stderr")
        mock_popen.return_value = mock_process

        inputs = [0.5] * 10
        flight_conditions = (0.5, 1000, 1.225, 50)
        result = self.loop_thread.run_datcom(inputs, flight_conditions, alpha_cruise=2.0, mass=1.1, prev_data_dict={}, to_trim=None)
        self.assertIsInstance(result, dict)

if __name__ == "__main__":
    unittest.main()