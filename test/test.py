# !/usr/bin/python
# -*- coding: utf-8 -*-
# Author:lihuiru
# Created on 2024/3/13 11:01
import argparse
import os
import subprocess
import sys
import traceback
import unittest
from contextlib import redirect_stdout, redirect_stderr

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class fluhpTest(unittest.TestCase):
    directory_path = "./"

    @classmethod
    def setUpClass(cls):
        # 在所有测试开始之前执行一次，用于设置测试所需的全局配置

        cls.test_files_path = os.path.join(cls.directory_path, 'test_files')
        cls.standardized_fasta_path = os.path.join(cls.directory_path, 'standardized_fasta')
        cls.adaptation_path = os.path.join(cls.directory_path, 'adaptation')

        cls.result_path = os.path.join(cls.directory_path, 'result')
        cls.ada_prediction_path = os.path.join(cls.directory_path, 'ada_prediction')


        # 创建所有需要的目录
        os.makedirs(cls.result_path, exist_ok = True)
        os.makedirs(cls.ada_prediction_path, exist_ok = True)
        # os.makedirs(cls.vir_prediction_path, exist_ok = True)

    def run_command(self, command):
        # 辅助函数，用于执行命令并返回结果 (modified to support Python 3.6)
        process = subprocess.run(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE,
                                 universal_newlines = True)
        return process

    def test_anno_command(self):
        # 测试anno功能
        command = f'fluhp anno -i {self.test_files_path} -o {self.result_path}'
        self.test_output = self.run_command(command)  # 修改这里，保存输出到self.test_output
        self.assertEqual(self.test_output.returncode, 0, "anno command failed with error: " + self.test_output.stderr)

    def test_extract_command(self):
        # 测试extract功能
        command = f'fluhp extract -i {self.standardized_fasta_path} -a {self.result_path}'
        self.test_output = self.run_command(command)
        self.assertEqual(self.test_output.returncode, 0,
                         "extract command failed with error: " + self.test_output.stderr)

    def test_predh_command(self):
        # 测试predh功能
        command = f'fluhp predh -i {self.adaptation_path} -o {self.ada_prediction_path}'
        self.test_output = self.run_command(command)
        self.assertEqual(self.test_output.returncode, 0, "predh command failed with error: " + self.test_output.stderr)


    def get_threshold(self):
        # 从命令行参数中获取阈值，默认为0.5
        parser = argparse.ArgumentParser()
        parser.add_argument('-th', '--threshold', default=0.5, type=float,
                            help='Probability threshold for model prediction.')
        args, _ = parser.parse_known_args()  # 解析命令行参数，忽略unittest无法识别的参数
        return args.threshold

    def test_t_create_detailed_info_zip(self):
        # 测试创建detailed_info.zip功能
        import zipfile
        import glob

        # 设置输出zip文件路径
        output_zip_path = os.path.join(self.directory_path, 'detailed_info.zip')

        # 查找所有以_markers.csv结尾的文件
        marker_files = glob.glob(os.path.join(self.adaptation_path, '*_markers.csv'))

        # 检查是否找到了文件
        self.assertTrue(len(marker_files) > 0, f"No *_markers.csv files found in {self.adaptation_path}")

        # 创建zip文件
        with zipfile.ZipFile(output_zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for file_path in marker_files:
                # 将文件添加到zip中，只使用文件名而不是完整路径
                zipf.write(file_path, os.path.basename(file_path))

        # 验证zip文件已创建并且不为空
        self.assertTrue(os.path.exists(output_zip_path), "detailed_info.zip was not created")
        self.assertTrue(os.path.getsize(output_zip_path) > 0, "detailed_info.zip is empty")

        print(f"Successfully created detailed_info.zip with {len(marker_files)} marker files")

        # 返回创建的zip文件路径，以便其他测试可能使用
        return output_zip_path

    def tearDown(self):
        # Method to perform cleanup after tests, including capturing logs and exceptions
        test_method_name = self.id().split('.')[-1]  # Getting the method name
        print(test_method_name)
        test_method_dir = os.path.join(self.directory_path, 'runlog', test_method_name)
        os.makedirs(test_method_dir, exist_ok = True)

        # If an exception occurred, capture the information
        exception_info = sys.exc_info()
        has_exception = exception_info[0] is not None

        log_path = os.path.join(test_method_dir, 'log.txt')

        # Capture the logs and exceptions (if any) in the log file
        with open(log_path, 'w') as log_file, redirect_stdout(log_file), redirect_stderr(log_file):
            if has_exception:
                traceback.print_exception(*exception_info, file = log_file)

            # Capture the stdout and stderr from the test method (if any)
            test_output = getattr(self, 'test_output', None)
            if test_output:
                stdout, stderr = test_output.stdout, test_output.stderr
                if stdout:
                    log_file.write("\nStandard Output:\n" + stdout)
                if stderr:
                    log_file.write("\nStandard Error:\n" + stderr)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-th', '--threshold', default=0.5, type=float,
                        help='Probability threshold for host prediction.')
    args, unittest_args = parser.parse_known_args()
    sys.argv[1:] = unittest_args
    unittest.main()
