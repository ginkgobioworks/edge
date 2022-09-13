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

    def test_no_duplicate_features(self):
        # Check that we reuse features and chunk features now
        feature1 = self.root.annotate(3, 5, "A1", "gene", 1, qualifiers={"test": 1})
        feature2 = self.root.annotate(3, 5, "A1", "gene", 1, qualifiers={"test": 1})
        self.assertEqual(feature1, feature2)
        self.assertEqual(
            list(feature1.chunk_feature_set.all()),
            list(feature2.chunk_feature_set.all())
        )

        # Check that we create new features and chunk features for different single values
        feature3 = self.root.annotate(3, 5, "A1", "gene", -1, qualifiers={"test": 1})
        self.assertNotEqual(feature1, feature3)
        self.assertNotEqual(
            list(feature1.chunk_feature_set.all()),
            list(feature3.chunk_feature_set.all())
        )
        self.assertEqual(
            feature1.chunk_feature_set.count(),
            feature3.chunk_feature_set.count()
        )

        # Check that we create new features for different qualifiers
        feature4 = self.root.annotate(3, 5, "A1", "gene", 1)
        self.assertNotEqual(feature1, feature4)

        # Check that we can force new feature generation
        feature5 = self.root.annotate(
            3, 5, "A1", "gene", 1,
            qualifiers={"test": 1}, force_new_annotation=True
        )
        self.assertNotEqual(feature1, feature5)
