from __future__ import annotations

import time
from typing import Any

import requests

from .models import Exon, SpeciesConfig, TranscriptRecord
from .sequence import clean_dna


SPECIES: dict[str, SpeciesConfig] = {
    "mouse": SpeciesConfig(
        key="mouse",
        label="Mouse (Mus musculus)",
        ensembl_name="mus_musculus",
        assembly="GRCm39",
    ),
    "human": SpeciesConfig(
        key="human",
        label="Human (Homo sapiens)",
        ensembl_name="homo_sapiens",
        assembly="GRCh38",
    ),
}


class EnsemblError(RuntimeError):
    pass


class EnsemblClient:
    """Small, explicit wrapper around the public Ensembl REST API."""

    def __init__(self, base_url: str = "https://rest.ensembl.org", timeout: float = 30.0):
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers.update(
            {
                "User-Agent": "HDRTagDesigner/0.3.1 (research prototype)",
                "Accept": "application/json",
            }
        )

    def _request_json(self, path: str, params: dict[str, Any] | None = None) -> Any:
        url = f"{self.base_url}{path}"
        try:
            response = self.session.get(
                url,
                params=params,
                headers={"Content-Type": "application/json"},
                timeout=self.timeout,
            )
        except requests.RequestException as exc:
            raise EnsemblError(f"Could not contact Ensembl REST: {exc}") from exc
        if response.status_code == 429:
            retry_after = float(response.headers.get("Retry-After", "1"))
            time.sleep(min(retry_after, 5.0))
            return self._request_json(path, params)
        if not response.ok:
            raise EnsemblError(
                f"Ensembl REST returned HTTP {response.status_code} for {path}: "
                f"{response.text[:300]}"
            )
        try:
            return response.json()
        except ValueError as exc:
            raise EnsemblError(f"Ensembl returned invalid JSON for {path}") from exc

    def _request_sequence(self, path: str, params: dict[str, Any] | None = None) -> str:
        url = f"{self.base_url}{path}"
        try:
            response = self.session.get(
                url,
                params=params,
                headers={"Content-Type": "text/plain", "Accept": "text/plain"},
                timeout=self.timeout,
            )
        except requests.RequestException as exc:
            raise EnsemblError(f"Could not contact Ensembl REST: {exc}") from exc
        if not response.ok:
            raise EnsemblError(
                f"Ensembl REST returned HTTP {response.status_code} for {path}: "
                f"{response.text[:300]}"
            )
        return clean_dna(response.text)

    def resolve_gene(self, species: SpeciesConfig, symbol_or_id: str) -> dict[str, Any]:
        query = symbol_or_id.strip()
        if query.upper().startswith(("ENSG", "ENSMUSG")):
            return self._request_json(f"/lookup/id/{query}", {"expand": 1})
        xrefs = self._request_json(
            f"/xrefs/symbol/{species.ensembl_name}/{query}",
            {"object_type": "gene"},
        )
        gene_hits = [item for item in xrefs if item.get("type") == "gene"]
        if not gene_hits:
            raise EnsemblError(f"No Ensembl gene found for {query!r} in {species.label}")
        # Prefer an exact display-label match, then use the first Ensembl result.
        exact = [
            item for item in gene_hits
            if str(item.get("display_id", "")).casefold() == query.casefold()
        ]
        gene_id = (exact or gene_hits)[0]["id"]
        return self._request_json(f"/lookup/id/{gene_id}", {"expand": 1})

    def transcript_record(
        self,
        species: SpeciesConfig,
        symbol_or_gene_id: str,
        transcript_id: str | None = None,
    ) -> TranscriptRecord:
        gene = self.resolve_gene(species, symbol_or_gene_id)
        transcripts = gene.get("Transcript", [])
        if not transcripts:
            raise EnsemblError(f"No transcripts returned for {gene.get('id', symbol_or_gene_id)}")

        requested = transcript_id.strip().split(".")[0] if transcript_id else ""
        canonical = str(gene.get("canonical_transcript", "")).split(".")[0]
        selected: dict[str, Any] | None = None
        if requested:
            selected = next((item for item in transcripts if item.get("id") == requested), None)
            if selected is None:
                # A direct transcript lookup handles transcripts not included in a reduced gene result.
                selected = self._request_json(f"/lookup/id/{requested}", {"expand": 1})
        elif canonical:
            selected = next((item for item in transcripts if item.get("id") == canonical), None)
        if selected is None:
            protein_coding = [
                item for item in transcripts if item.get("biotype") == "protein_coding"
            ]
            selected = (protein_coding or transcripts)[0]

        stable_id = selected["id"]
        expanded = self._request_json(f"/lookup/id/{stable_id}", {"expand": 1})
        cdna = self._request_sequence(f"/sequence/id/{stable_id}", {"type": "cdna"})
        cds = self._request_sequence(f"/sequence/id/{stable_id}", {"type": "cds"})
        exons = [
            Exon(
                start1=int(exon["start"]),
                end1=int(exon["end"]),
                strand=int(exon.get("strand", expanded.get("strand", 1))),
                stable_id=str(exon.get("id", "")),
            )
            for exon in expanded.get("Exon", [])
        ]
        if not exons:
            raise EnsemblError(f"No exon coordinates returned for {stable_id}")

        return TranscriptRecord(
            species=species,
            gene_symbol=str(gene.get("display_name") or symbol_or_gene_id),
            gene_id=str(gene["id"]),
            transcript_id=stable_id,
            transcript_version=str(expanded.get("version", "")),
            chromosome=str(expanded.get("seq_region_name", gene.get("seq_region_name", ""))),
            strand=int(expanded.get("strand", gene.get("strand", 1))),
            cdna=cdna,
            cds=cds,
            exons=exons,
            source="Ensembl REST",
            source_release="live",
        )

    def region_sequence(
        self,
        species: SpeciesConfig,
        chromosome: str,
        start0: int,
        end0: int,
    ) -> str:
        if start0 < 0 or end0 <= start0:
            raise ValueError(f"Invalid genomic interval [{start0}, {end0})")
        # Ensembl region endpoint uses 1-based inclusive coordinates.
        region = f"{chromosome}:{start0 + 1}..{end0}:1"
        return self._request_sequence(
            f"/sequence/region/{species.ensembl_name}/{region}",
        )
