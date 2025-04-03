import unittest
from unittest.mock import patch, MagicMock
from scripts.gui import TkinterInterface
import tkinter as tk

class TestTkinterInterface(unittest.TestCase):
    @patch("scripts.gui.AerodynamicOptimization")
    def test_run_results(self, mock_optimizer):
        root = tk.Tk()
        app = TkinterInterface(root)

        # Mock AerodynamicOptimization
        mock_instance = MagicMock()
        mock_optimizer.return_value = mock_instance

        # Simulate running results
        app.run_results()

        # Check if the optimization thread was started
        self.assertTrue(mock_optimizer.called)

if __name__ == "__main__":
    unittest.main()