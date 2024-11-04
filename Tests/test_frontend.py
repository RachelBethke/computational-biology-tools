import unittest
import tkinter as tk
import os
from src.genetics.frontend import (
    load_haplotype_data,
    process_block,
    handle_file_error,
    load_file,
    create_plot,
)

class TestFrontend(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls): #root tesing window
        cls.root = tk.Tk()
        cls.root.withdraw()  #hide the window (?)
        cls.sample_blocks = {
            0: [1, 2, 3],
            1: [4, 5, 6],
            2: [7, 8, 9],
            3: [10, 11, 12]
        }
        cls.sample_file_path = "data/haplotype/frontend_testing_data.txt"
        os.makedirs(os.path.dirname(cls.sample_file_path), exist_ok=True)


    def test_load_haplotype_data(self):
        # make fake file path
        with open(self.sample_file_path, "w") as file:
            for i, block in enumerate(self.sample_blocks.values()):
                file.write(",".join(map(str, block)))
                if i < len(self.sample_blocks) - 1:
                    file.write("\n\n")
        # with open("frontend_testing_data.txt", "w") as file:
        #     #file.write("\n\n".join(self.sample_blocks))
        #     file.write("\n\n".join(str(v) for v in self.sample_blocks))
        
        blocks = load_haplotype_data(self.sample_file_path)
        self.assertEqual(len(blocks), 4) 
        self.assertEqual(blocks[0].strip(), "1,2,3")
    
    def test_process_block(self):
        """
        Test processing a valid block of haplotype data
        """
        block = "1,2,3\n4,5,6"
        df = process_block(block)
        self.assertEqual(df.shape, (2, 3))
        self.assertTrue((df.values == [[1, 2, 3], [4, 5, 6]]).all())

    def test_empty_block(self):
        """
        Test the processing of an empty block
        """
        with self.assertRaises(ValueError):
            process_block("\n\n")
    
    def test_file_errors(self):
        """
        Test error handling
        """
        error_message = None
        try:
            raise ZeroDivisionError("division by zero")
        except ZeroDivisionError as e:
            error_message = str(e)
            handle_file_error(e)
        self.assertIn("division by zero", error_message)
    
    def test_elements(self):
        """
        Test main UI elements being created
        """
        main_frame = tk.Frame(self.root)
        upload_button = tk.Button(main_frame, text="Upload Haplo File", command=lambda: None)
        result_label = tk.Label(main_frame, text="Results here")
        
        self.assertEqual(upload_button.cget("text"), "Upload Haplo File")
        self.assertEqual(result_label.cget("text"), "Results here")
    
    @classmethod
    def tearDownClass(cls):
        cls.root.destroy()

if __name__ == "__main__":
    unittest.main()
