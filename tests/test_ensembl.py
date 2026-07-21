from __future__ import annotations

import unittest
from unittest.mock import patch

import requests

from hdr_designer.ensembl import EnsemblClient, EnsemblError, SPECIES


class FakeResponse:
    def __init__(
        self,
        status_code: int,
        *,
        json_value=None,
        text: str = "",
        headers: dict[str, str] | None = None,
    ) -> None:
        self.status_code = status_code
        self._json_value = json_value
        self.text = text
        self.headers = headers or {}

    @property
    def ok(self) -> bool:
        return 200 <= self.status_code < 400

    def json(self):
        return self._json_value


class FakeSession:
    def __init__(self, responses: list[FakeResponse | Exception]) -> None:
        self.responses = list(responses)
        self.calls = 0
        self.headers: dict[str, str] = {}

    def get(self, *args, **kwargs):
        self.calls += 1
        response = self.responses.pop(0)
        if isinstance(response, Exception):
            raise response
        return response


class EnsemblRetryTest(unittest.TestCase):
    @patch("hdr_designer.ensembl.time.sleep")
    def test_429_retries_then_returns_json(self, sleep) -> None:
        client = EnsemblClient(max_retries=2, retry_backoff=0)
        session = FakeSession(
            [
                FakeResponse(429, headers={"Retry-After": "0"}),
                FakeResponse(200, json_value={"id": "ENSG1"}),
            ]
        )
        client.session = session
        self.assertEqual(client._request_json("/lookup/id/ENSG1"), {"id": "ENSG1"})
        self.assertEqual(session.calls, 2)
        sleep.assert_called_once_with(0.0)

    @patch("hdr_designer.ensembl.time.sleep")
    def test_500_retry_is_bounded_and_error_omits_html(self, sleep) -> None:
        client = EnsemblClient(max_retries=2, retry_backoff=0)
        session = FakeSession(
            [
                FakeResponse(500, text="<html>server details</html>"),
                FakeResponse(503, text="<html>server details</html>"),
                FakeResponse(500, text="<html>server details</html>"),
            ]
        )
        client.session = session
        with self.assertRaisesRegex(EnsemblError, "temporarily unavailable") as raised:
            client._request_json("/lookup/id/ENST1")
        self.assertEqual(session.calls, 3)
        self.assertNotIn("<html>", str(raised.exception))
        self.assertEqual(sleep.call_count, 2)

    @patch("hdr_designer.ensembl.time.sleep")
    def test_sequence_request_retries_network_failure(self, sleep) -> None:
        client = EnsemblClient(max_retries=1, retry_backoff=0)
        session = FakeSession(
            [
                requests.ConnectionError("temporary"),
                FakeResponse(200, text="ACGT\n"),
            ]
        )
        client.session = session
        self.assertEqual(client._request_sequence("/sequence/id/ENST1"), "ACGT")
        self.assertEqual(session.calls, 2)

    def test_supplied_transcript_must_belong_to_selected_gene(self) -> None:
        client = EnsemblClient(max_retries=0)
        gene = {
            "id": "ENSG_SELECTED",
            "display_name": "GENE",
            "Transcript": [{"id": "ENST_OTHER", "biotype": "protein_coding"}],
        }
        with patch.object(client, "resolve_gene", return_value=gene), patch.object(
            client,
            "_request_json",
            return_value={
                "id": "ENST_REQUESTED",
                "Parent": "ENSG_DIFFERENT",
                "biotype": "protein_coding",
            },
        ):
            with self.assertRaisesRegex(EnsemblError, "belongs to gene ENSG_DIFFERENT"):
                client.transcript_record(
                    SPECIES["human"],
                    "GENE",
                    "ENST_REQUESTED",
                )


if __name__ == "__main__":
    unittest.main()
