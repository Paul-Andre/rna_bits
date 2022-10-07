import unittest

import parser


class ParseParensTestCase(unittest.TestCase):
    def test_success(self):
        cases = [
            ("", set()),
            (".", set()),
            ("..", set()),
            ("()", {(1, 2)}),
            ("{}", {(1, 2)}),
            ("<>", {(1, 2)}),
            ("[]", {(1, 2)}),
            ("Aa", {(1, 2)}),
            ("Bb", {(1, 2)}),
            ("([)]", {(1, 3), (2, 4)}),
            ("(&)", {(1, 2)}),
            ("(.)", {(1, 3)}),
            ("(.())", {(1, 5), (3, 4)}),
        ]
        for input, expected_output in cases:
            actual_output = set(parser.parseParens(input))
            self.assertEqual(expected_output, actual_output, "For test " + repr(input))


class SplitAndPositionsTestCase(unittest.TestCase):
    def test_success(self):
        self.assertEqual(
            list(parser.splitAndPositions("asdf: sd", ":")), [(0, "asdf"), (5, " sd")]
        )
        self.assertEqual(
            list(parser.splitAndPositions("asdf: sd", ":", 100)),
            [(100, "asdf"), (105, " sd")],
        )


if __name__ == "__main__":
    unittest.main()
