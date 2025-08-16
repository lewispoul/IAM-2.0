# Ketcher PostMessage Contract

**Sources:** IAM codebase, Ketcher HTML/JS bridge  
**Note:** Only the first 10 search results are shown. [View more search results in GitHub](https://github.com/lewispoul/IAM/search?q=postMessage+ketcher+iframe+message+mol+smiles+origin&type=code)

---

## Contract Overview

Communication between the parent IAM app and the embedded Ketcher iframe (structure editor) uses the browser's `window.postMessage` API.  
**Security:** Messages are strictly filtered to/from the same origin (`event.origin === window.location.origin`).  
**Payloads:** All messages use plain JSON objects.

---

## Message Types Sent **From Parent to Ketcher iframe**

### 1. `get-molfile`
- **Purpose:** Request Ketcher to send the current molecule in MOL format.
- **Schema:**
  ```json
  {
    "type": "get-molfile"
  }
  ```
- **Example:**
  ```json
  { "type": "get-molfile" }
  ```
- **Sent via:**  
  ```js
  ketcherFrame.contentWindow.postMessage({ type: 'get-molfile' }, '*');
  ```

### 2. (Potential future) `set-molecule`
- **Purpose:** Set the molecule in Ketcher (not shown in main code, but typical in sketcher contracts).
- **Schema:**
  ```json
  {
    "type": "set-molecule",
    "data": "<molfile or SMILES string>",
    "format": "mol" | "smiles"
  }
  ```
- **Example:**
  ```json
  { "type": "set-molecule", "data": "C1=CC=CC=C1", "format": "smiles" }
  ```

---

## Message Types Sent **From Ketcher iframe to Parent**

### 1. `molfile`
- **Purpose:** Send the current molecule's MOL file to the parent.
- **Schema:**
  ```json
  {
    "type": "molfile",
    "molfile": "<mol block as string>"
  }
  ```
- **Example:**
  ```json
  { "type": "molfile", "molfile": "\n  Ketcher 3.4.0 ...\n..." }
  ```
- **Sent via:**  
  ```js
  window.parent.postMessage({ type: 'molfile', molfile }, event.origin);
  ```

### 2. (Potential future) `smiles`
- **Purpose:** Send the current molecule as SMILES to the parent.
- **Schema:**
  ```json
  {
    "type": "smiles",
    "smiles": "<SMILES string>"
  }
  ```
- **Example:**
  ```json
  { "type": "smiles", "smiles": "C1=CC=CC=C1" }
  ```

---

## Security Checks

- **Origin Filtering:**  
  Both parent and iframe listeners check:
  ```js
  if (event.origin !== window.location.origin) return;
  ```
  Only messages from the exact same origin are processed, avoiding cross-site scripting attacks.

- **Source Filtering (parent side):**
  ```js
  if (!event.data || event.source !== ketcherFrame.contentWindow) return;
  ```
  Ensures only the intended iframe responds.

---

## Example Message Exchange

**Parent → Ketcher:**
```json
{ "type": "get-molfile" }
```

**Ketcher → Parent:**
```json
{ "type": "molfile", "molfile": "Mol block content..." }
```

---

## Full JSON Schemas

```json
{
  "get-molfile": {
    "type": "object",
    "properties": {
      "type": { "const": "get-molfile" }
    },
    "required": ["type"]
  },
  "molfile": {
    "type": "object",
    "properties": {
      "type": { "const": "molfile" },
      "molfile": { "type": "string" }
    },
    "required": ["type", "molfile"]
  },
  "set-molecule": {
    "type": "object",
    "properties": {
      "type": { "const": "set-molecule" },
      "data": { "type": "string" },
      "format": { "enum": ["mol", "smiles"] }
    },
    "required": ["type", "data", "format"]
  },
  "smiles": {
    "type": "object",
    "properties": {
      "type": { "const": "smiles" },
      "smiles": { "type": "string" }
    },
    "required": ["type", "smiles"]
  }
}
```

---

## Markdown Spec

### IAM–Ketcher iframe PostMessage Contract

- All messages are JSON objects and must include a `type` field.
- **Parent → iframe:**
  - `get-molfile`: Requests MOL block for current molecule.
  - (optional) `set-molecule`: Sets molecule in Ketcher by MOL or SMILES.
- **iframe → Parent:**
  - `molfile`: Sends MOL block for current molecule.
  - (optional) `smiles`: Sends SMILES for current molecule.
- **Security:** Only process messages where `event.origin === window.location.origin` and correct `event.source`.

---

**For more integration details, see**  
[GitHub code search results](https://github.com/lewispoul/IAM/search?q=postMessage+ketcher+iframe+message+mol+smiles+origin&type=code)
