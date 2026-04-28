/**
 * Data Downloads page — manage vendor sessions, kick off downloads,
 * watch progress.
 *
 * MVP: polls /jobs every 5s. WebSocket realtime push is wired up in
 * the backend (realtime:task:<job_id>); a follow-up can subscribe via
 * the existing useRealtimeConnection hook to remove the polling.
 */
import { useMemo, useState } from 'react'
import {
  Button,
  Card,
  Form,
  Input,
  Modal,
  Progress,
  Select,
  Space,
  Table,
  Tag,
  Tooltip,
  Typography,
  message,
  Switch,
} from 'antd'
import { DeleteOutlined, ReloadOutlined, KeyOutlined, PlayCircleOutlined } from '@ant-design/icons'
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query'
import { useDownloadJobsRealtime } from '@/hooks/useRealtime'
import { dataDownloadService, vendorCredentialService } from '@/services/data-download.service'
import type { DownloadJob, DownloadJobCreateRequest, DownloadStatus, SessionLoginRequest } from '@/types/data-download'

const { Title, Text, Paragraph } = Typography

const STATUS_COLOR: Record<DownloadStatus, string> = {
  pending: 'default',
  running: 'processing',
  completed: 'success',
  failed: 'error',
  cancelled: 'warning',
}

const STATUS_LABEL: Record<DownloadStatus, string> = {
  pending: '等待中',
  running: '下载中',
  completed: '已完成',
  failed: '失败',
  cancelled: '已取消',
}

// Vendor id → display name. Falls back to id for unknown vendors.
const VENDOR_LABEL: Record<string, string> = {
  lc_bio: '联川生物 (lcbio)',
  novogene: '诺禾致源 (即将支持)',
}

const formatVendor = (id: string) => VENDOR_LABEL[id] ?? id

export function DataDownloadsPage() {
  const qc = useQueryClient()
  const [vendor, setVendor] = useState<string>('lc_bio')
  const [credentialModalOpen, setCredentialModalOpen] = useState(false)

  // Catalog of vendors the backend has adapters for.
  const vendorsQuery = useQuery({
    queryKey: ['data-downloads', 'vendors'],
    queryFn: () => dataDownloadService.listVendors(),
  })

  // Saved credentials for the chosen vendor.
  const credsQuery = useQuery({
    queryKey: ['vendor-credentials', vendor],
    queryFn: () => vendorCredentialService.list(vendor).then((r) => r.items),
  })

  // Whether the vendor's daemon is currently up.
  const sessionQuery = useQuery({
    queryKey: ['data-downloads', 'session', vendor],
    queryFn: () => dataDownloadService.getSession(vendor),
    refetchInterval: 10_000,
  })

  // User's download jobs. Realtime push (useDownloadJobsRealtime below)
  // invalidates the cache on each progress event; the polling interval
  // is just a slow safety net for missed events.
  const jobsQuery = useQuery({
    queryKey: ['data-downloads', 'jobs'],
    queryFn: () => dataDownloadService.listJobs().then((r) => r.items),
    refetchInterval: (query) => {
      const items = (query.state.data ?? []) as DownloadJob[]
      return items.some((j) => j.status === 'running' || j.status === 'pending') ? 30_000 : false
    },
  })

  // Subscribe to realtime progress for any in-flight job so the table
  // updates within ~100ms of each lcbio progress line.
  const inFlightIds = useMemo(
    () => (jobsQuery.data ?? []).filter((j) => j.status === 'running' || j.status === 'pending').map((j) => j.id),
    [jobsQuery.data],
  )
  useDownloadJobsRealtime(inFlightIds)

  const openSessionMutation = useMutation({
    mutationFn: (req: SessionLoginRequest) => dataDownloadService.openSession(req),
    onSuccess: () => {
      message.success('会话已开启')
      qc.invalidateQueries({ queryKey: ['data-downloads', 'session', vendor] })
    },
    onError: (e: Error & { response?: { data?: { detail?: string } } }) =>
      message.error(e.response?.data?.detail ?? e.message),
  })

  const closeSessionMutation = useMutation({
    mutationFn: () => dataDownloadService.closeSession(vendor),
    onSuccess: () => {
      message.success('会话已关闭')
      qc.invalidateQueries({ queryKey: ['data-downloads', 'session', vendor] })
    },
  })

  const createJobMutation = useMutation({
    mutationFn: (req: DownloadJobCreateRequest) => dataDownloadService.createJob(req),
    onSuccess: () => {
      message.success('下载任务已创建')
      qc.invalidateQueries({ queryKey: ['data-downloads', 'jobs'] })
    },
    onError: (e: Error & { response?: { data?: { detail?: string } } }) =>
      message.error(e.response?.data?.detail ?? e.message),
  })

  const cancelJobMutation = useMutation({
    mutationFn: (id: string) => dataDownloadService.cancelJob(id),
    onSuccess: () => qc.invalidateQueries({ queryKey: ['data-downloads', 'jobs'] }),
  })

  const sessionActive = sessionQuery.data?.active === true

  return (
    <div style={{ padding: 24 }}>
      <Title level={3}>数据下载</Title>
      <Paragraph type="secondary">
        从测序公司（联川 / 诺禾致源 / ...）直接将原始数据下载到 NGSmodule。 在公司网页的下载按钮上复制 obs 路径（
        <Text code>download</Text> 后引号里的部分），粘贴到下面即可。
      </Paragraph>

      <Card
        title={
          <Space>
            <span>厂商</span>
            <Select
              size="small"
              style={{ minWidth: 200 }}
              loading={vendorsQuery.isLoading}
              value={vendor}
              onChange={setVendor}
              options={(vendorsQuery.data ?? []).map((v: string) => ({ value: v, label: formatVendor(v) }))}
            />
            <Tag color={sessionActive ? 'success' : 'default'}>{sessionActive ? '会话已开启' : '会话未开启'}</Tag>
          </Space>
        }
        extra={
          <Space>
            <Button icon={<KeyOutlined />} onClick={() => setCredentialModalOpen(true)}>
              已保存登录（{credsQuery.data?.length ?? 0}）
            </Button>
            {sessionActive && (
              <Button danger onClick={() => closeSessionMutation.mutate()}>
                关闭会话
              </Button>
            )}
          </Space>
        }
        style={{ marginBottom: 16 }}
      >
        <SessionForm
          vendor={vendor}
          credentials={credsQuery.data ?? []}
          onSubmit={(req) => openSessionMutation.mutate(req)}
          submitting={openSessionMutation.isPending}
          disabled={sessionActive}
        />
      </Card>

      <Card title="新建下载" style={{ marginBottom: 16 }}>
        <DownloadForm
          vendor={vendor}
          disabled={!sessionActive}
          submitting={createJobMutation.isPending}
          onSubmit={(req) => createJobMutation.mutate(req)}
        />
      </Card>

      <Card
        title={`下载任务（${jobsQuery.data?.length ?? 0}）`}
        extra={
          <Button
            icon={<ReloadOutlined />}
            onClick={() => qc.invalidateQueries({ queryKey: ['data-downloads', 'jobs'] })}
          >
            刷新
          </Button>
        }
      >
        <JobsTable
          jobs={jobsQuery.data ?? []}
          loading={jobsQuery.isLoading}
          onCancel={(id) => cancelJobMutation.mutate(id)}
        />
      </Card>

      <CredentialsModal vendor={vendor} open={credentialModalOpen} onClose={() => setCredentialModalOpen(false)} />
    </div>
  )
}

// -------- session form --------

function SessionForm({
  vendor,
  credentials,
  onSubmit,
  submitting,
  disabled,
}: {
  vendor: string
  credentials: { id: string; label: string; email_preview: string }[]
  onSubmit: (req: SessionLoginRequest) => void
  submitting: boolean
  disabled: boolean
}) {
  const [credentialId, setCredentialId] = useState<string | undefined>(undefined)
  const [form] = Form.useForm<{ email: string; password: string }>()

  return (
    <Form
      layout="inline"
      form={form}
      disabled={disabled}
      onFinish={(values) => {
        if (credentialId) {
          onSubmit({ vendor, credential_id: credentialId })
        } else {
          onSubmit({ vendor, ...values })
        }
      }}
    >
      <Form.Item label="使用已保存">
        <Select
          allowClear
          style={{ minWidth: 220 }}
          placeholder="（无 — 下方手动输入）"
          value={credentialId}
          onChange={setCredentialId}
          options={credentials.map((c) => ({
            value: c.id,
            label: `${c.label} — ${c.email_preview}`,
          }))}
        />
      </Form.Item>
      {!credentialId && (
        <>
          <Form.Item name="email" rules={[{ required: !credentialId, message: '请填写邮箱' }]}>
            <Input placeholder="邮箱" autoComplete="off" style={{ width: 220 }} />
          </Form.Item>
          <Form.Item name="password" rules={[{ required: !credentialId, message: '请填写密码' }]}>
            <Input.Password placeholder="密码" autoComplete="off" style={{ width: 200 }} />
          </Form.Item>
        </>
      )}
      <Form.Item>
        <Button type="primary" htmlType="submit" icon={<PlayCircleOutlined />} loading={submitting}>
          开启会话
        </Button>
      </Form.Item>
    </Form>
  )
}

// -------- new download form --------

function DownloadForm({
  vendor,
  disabled,
  submitting,
  onSubmit,
}: {
  vendor: string
  disabled: boolean
  submitting: boolean
  onSubmit: (req: DownloadJobCreateRequest) => void
}) {
  const [form] = Form.useForm<{
    source_path: string
    dest_path: string
    auto_register: boolean
    project_name?: string
  }>()
  return (
    <Form
      form={form}
      layout="vertical"
      disabled={disabled}
      onFinish={(values) =>
        onSubmit({
          vendor,
          source_path: values.source_path,
          dest_path: values.dest_path,
          auto_register: values.auto_register,
          project_name: values.project_name,
        })
      }
      initialValues={{ auto_register: false }}
    >
      <Form.Item
        name="source_path"
        label="源路径（公司网页复制的 obs 路径）"
        rules={[{ required: true, message: '请粘贴 obs 路径' }]}
      >
        <Input placeholder="LC-P20260327047-.../2026-04-21_17:24:55/Data.tar" />
      </Form.Item>
      <Form.Item name="dest_path" label="本地保存目录" rules={[{ required: true, message: '请填写保存路径' }]}>
        <Input placeholder="/data/lab/.../rawdata" />
      </Form.Item>
      <Form.Item name="auto_register" label="完成后自动注册为 NGSmodule 项目" valuePropName="checked">
        <Switch />
      </Form.Item>
      <Form.Item name="project_name" label="项目名称（可选，留空则按源路径自动推断）">
        <Input maxLength={100} />
      </Form.Item>
      <Form.Item>
        <Button type="primary" htmlType="submit" loading={submitting} disabled={disabled}>
          开始下载
        </Button>
      </Form.Item>
    </Form>
  )
}

// -------- jobs table --------

function JobsTable({
  jobs,
  loading,
  onCancel,
}: {
  jobs: DownloadJob[]
  loading: boolean
  onCancel: (id: string) => void
}) {
  const columns = useMemo(
    () => [
      {
        title: '源文件',
        dataIndex: 'source_path',
        key: 'source_path',
        ellipsis: true,
        render: (s: string) => (
          <Tooltip title={s}>
            <Text code>{s.split('/').pop()}</Text>
          </Tooltip>
        ),
      },
      {
        title: '保存目录',
        dataIndex: 'dest_path',
        key: 'dest_path',
        ellipsis: true,
        width: 220,
      },
      {
        title: '状态',
        dataIndex: 'status',
        key: 'status',
        width: 110,
        render: (s: DownloadStatus) => <Tag color={STATUS_COLOR[s]}>{STATUS_LABEL[s]}</Tag>,
      },
      {
        title: '进度',
        dataIndex: 'progress_pct',
        key: 'progress_pct',
        width: 200,
        render: (pct: number, row: DownloadJob) => (
          <Progress
            percent={Math.round(pct)}
            status={row.status === 'failed' ? 'exception' : row.status === 'completed' ? 'success' : 'active'}
            size="small"
          />
        ),
      },
      {
        title: '开始时间',
        dataIndex: 'started_at',
        key: 'started_at',
        width: 170,
        render: (v: string | null) => (v ? new Date(v).toLocaleString() : '—'),
      },
      {
        title: '',
        key: 'actions',
        width: 80,
        render: (_: unknown, row: DownloadJob) =>
          row.status === 'running' || row.status === 'pending' ? (
            <Button size="small" danger icon={<DeleteOutlined />} onClick={() => onCancel(row.id)}>
              取消
            </Button>
          ) : null,
      },
    ],
    [onCancel],
  )
  return (
    <Table
      rowKey="id"
      size="small"
      loading={loading}
      dataSource={jobs}
      columns={columns}
      pagination={{ pageSize: 10 }}
    />
  )
}

// -------- saved credentials modal --------

function CredentialsModal({ vendor, open, onClose }: { vendor: string; open: boolean; onClose: () => void }) {
  const qc = useQueryClient()
  const [form] = Form.useForm<{ email: string; password: string; label: string }>()

  const credsQuery = useQuery({
    queryKey: ['vendor-credentials', vendor],
    queryFn: () => vendorCredentialService.list(vendor).then((r) => r.items),
    enabled: open,
  })

  const createMutation = useMutation({
    mutationFn: (req: { vendor: string; email: string; password: string; label: string }) =>
      vendorCredentialService.create(req),
    onSuccess: () => {
      message.success('已保存')
      form.resetFields()
      qc.invalidateQueries({ queryKey: ['vendor-credentials', vendor] })
    },
  })

  const deleteMutation = useMutation({
    mutationFn: (id: string) => vendorCredentialService.delete(id),
    onSuccess: () => qc.invalidateQueries({ queryKey: ['vendor-credentials', vendor] }),
  })

  return (
    <Modal
      open={open}
      onCancel={onClose}
      title={`已保存的登录信息（${formatVendor(vendor)}）`}
      footer={null}
      width={640}
    >
      <Table
        rowKey="id"
        size="small"
        dataSource={credsQuery.data ?? []}
        pagination={false}
        locale={{ emptyText: '暂无保存的登录信息' }}
        columns={[
          { title: '标签', dataIndex: 'label' },
          { title: '邮箱', dataIndex: 'email_preview' },
          {
            title: '',
            width: 70,
            render: (_, row) => (
              <Button size="small" danger icon={<DeleteOutlined />} onClick={() => deleteMutation.mutate(row.id)} />
            ),
          },
        ]}
        style={{ marginBottom: 16 }}
      />
      <Form
        form={form}
        layout="inline"
        onFinish={(v) => createMutation.mutate({ vendor, ...v, label: v.label || 'default' })}
      >
        <Form.Item name="label">
          <Input placeholder="标签（默认 default）" style={{ width: 150 }} />
        </Form.Item>
        <Form.Item name="email" rules={[{ required: true, message: '邮箱必填' }]}>
          <Input placeholder="邮箱" style={{ width: 180 }} autoComplete="off" />
        </Form.Item>
        <Form.Item name="password" rules={[{ required: true, message: '密码必填' }]}>
          <Input.Password placeholder="密码" style={{ width: 160 }} autoComplete="off" />
        </Form.Item>
        <Form.Item>
          <Button type="primary" htmlType="submit" loading={createMutation.isPending}>
            保存
          </Button>
        </Form.Item>
      </Form>
    </Modal>
  )
}
