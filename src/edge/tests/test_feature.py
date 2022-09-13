from django.test import TestCase

from edge.models import Fragment


class FeatureTests(TestCase):
    def setUp(self):
        self.root_sequence = "agttcgaggctga"
        self.root = Fragment.create_with_sequence("Foo", self.root_sequence)

    def test_sequence_positive_strand(self):
        feature = self.root.annotate(3, 5, "A1", "gene", 1)
        self.assertEqual(feature.sequence, "ttc")

    def test_sequence_negative_strand(self):
        feature = self.root.annotate(3, 5, "A1", "gene", -1)
        self.assertEqual(feature.sequence, "gaa")
