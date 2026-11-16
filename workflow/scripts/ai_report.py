"""Optional AI step: Claude drafts the Results section and figure legends from key_numbers.json.
The model only sees numbers the pipeline computed, is told to cite nothing else, and the draft is marked for expert review.
Credentials: ANTHROPIC_API_KEY, or an `ant auth login` profile (the SDK resolves either). --dry-run prints the prompt only."""
import argparse, json, sys

ap = argparse.ArgumentParser()
ap.add_argument("--numbers", required=True); ap.add_argument("--config", required=True, help="JSON: species, focal, tools")
ap.add_argument("--model", default="claude-opus-5-5"); ap.add_argument("--effort", default="high")
ap.add_argument("--out", required=True); ap.add_argument("--dry-run", action="store_true")
A = ap.parse_args()
numbers = json.load(open(A.numbers)); cfg = json.loads(A.config)

SYSTEM = """You are drafting the Results section of a comparative-genomics paper on lactic acid bacteria for a food-microbiology group.
Write for expert readers in the style of BMC Genomics: past tense, precise, no speculation beyond what the numbers support.
Use only the numbers in the provided JSON. Do not invent values, strains, genes, p-values or references. If a section has no data, omit it.
When an interpretation goes beyond the data (for example, an ecological explanation), label it as a hypothesis.
Output Markdown: one '##' subsection per analysis (dataset and species structure; pangenome; phylogeny and the focal strain;
focal-strain gene gain and loss; bacteriocins; prophages; CRISPR-Cas; CAZymes; plasmids; insertion sequences; safety),
then '## Figure legends' with one short legend per figure named in the JSON 'figures' list."""

user = (f"Species: {cfg['species']}. Focal strain: {numbers['focal']['strain']} ({numbers['focal']['accession']}), isolated from {numbers['focal']['niche']}.\n"
        f"Tools run (group pipeline): {', '.join(cfg['tools'])}.\nFigures produced: {', '.join(cfg['figures'])}.\n\n"
        f"Computed results (JSON):\n```json\n{json.dumps(numbers, indent=1, default=str)}\n```")
if A.dry_run:
    print(SYSTEM, "\n\n", user); sys.exit(0)

import anthropic  # imported here so --dry-run works without the SDK

client = anthropic.Anthropic()
try:
    # server-side fallback re-runs a refused request on another model inside the same call
    with client.beta.messages.stream(
        model=A.model, max_tokens=64000, system=SYSTEM,
        messages=[{"role": "user", "content": user}],
        output_config={"effort": A.effort},
        betas=["server-side-fallback-2026-07-01"], fallbacks="default",
    ) as stream:
        msg = stream.get_final_message()
except anthropic.AuthenticationError:
    sys.exit("No valid Anthropic credentials: set ANTHROPIC_API_KEY or run `ant auth login`.")
except anthropic.RateLimitError:
    sys.exit("Rate limited by the Anthropic API; rerun later.")
except anthropic.APIStatusError as e:
    sys.exit(f"Anthropic API error {e.status_code} (request {e.response.headers.get('request-id')}).")
except anthropic.APIConnectionError:
    sys.exit("Could not reach the Anthropic API.")
if msg.stop_reason == "refusal":
    sys.exit("The model declined to write the report; see stop_details in the logs.")
text = "".join(b.text for b in msg.content if b.type == "text")
note = ("> DRAFT written by Claude from the pipeline's computed numbers (results/report/key_numbers.json). "
        "Check every statement against the tables before use.\n\n")
open(A.out, "w").write(note + text + ("\n\n_[truncated at max_tokens]_\n" if msg.stop_reason == "max_tokens" else ""))
print(f"draft written to {A.out} ({msg.usage.output_tokens} output tokens, model {msg.model})")
