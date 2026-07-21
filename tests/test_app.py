from __future__ import annotations

from pathlib import Path
import unittest

from streamlit.testing.v1 import AppTest


class StreamlitAppTest(unittest.TestCase):
    def test_default_tubb5_workflow(self) -> None:
        app_path = Path(__file__).resolve().parents[1] / "app.py"
        app = AppTest.from_file(str(app_path), default_timeout=20).run()
        self.assertEqual(len(app.exception), 0)
        self.assertEqual(len(app.button), 1)

        app.button[0].click().run()
        self.assertEqual(len(app.exception), 0)
        self.assertEqual(
            [message.value for message in app.success],
            [
                "SEQUENCE-COMPLETE COMPUTATIONAL DESIGN",
                "All 1 internal SapI site(s) were removed from the final arm sequences.",
                "All three genotyping assays have a primer pair.",
            ],
        )
        metrics = {metric.label: metric.value for metric in app.metric}
        self.assertEqual(metrics["Transcript"], "ENSMUST00000001566.10")
        self.assertEqual(metrics["Fusion"], "686 aa")
        self.assertEqual(metrics["Final plasmid"], "3950 bp")
        self.assertEqual(metrics["Final SapI sites"], "0")
        self.assertEqual(metrics["Arm SapI sites found"], "1")
        self.assertEqual(metrics["Arm SapI sites resolved"], "1")
        self.assertEqual(metrics["Arm SapI sites remaining"], "0")
        self.assertEqual(len(app.get("download_button")), 8)

        # A normal Streamlit rerun (including older download behavior) must
        # restore the completed result from session state instead of resetting.
        app.run()
        self.assertEqual(len(app.exception), 0)
        self.assertIn(
            "SEQUENCE-COMPLETE COMPUTATIONAL DESIGN",
            [message.value for message in app.success],
        )
        self.assertEqual(len(app.get("download_button")), 8)

    def test_n_terminal_selection_uses_169226_and_disables_fixture(self) -> None:
        app_path = Path(__file__).resolve().parents[1] / "app.py"
        app = AppTest.from_file(str(app_path), default_timeout=20).run()
        app.selectbox[1].select("N-terminal").run()
        self.assertEqual(len(app.exception), 0)
        text_inputs = {item.label: item for item in app.text_input}
        self.assertEqual(
            text_inputs["Selected backbone"].value,
            "TVBB N-term-mNeongreen (Addgene #169226)",
        )
        fixture = next(
            item
            for item in app.checkbox
            if item.label == "Use bundled Tubb5 validation fixture"
        )
        self.assertFalse(fixture.value)
        self.assertTrue(fixture.disabled)
        self.assertFalse(text_inputs["3. Gene symbol or Ensembl gene ID"].disabled)


if __name__ == "__main__":
    unittest.main()
